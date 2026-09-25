# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: animate_gpp_robinson.R
# Description: Generates the Gross Primary Productivity climatology
#     animation from MOD17A2HGF (2000-2009), rendered in the Robinson
#     projection.
# Note: For the Plate Carree, Foucaut, and Nicolosi projections, see the other
#     scripts in the gpp scripts folder. Data excludes 2000-01-01 through
#     2000-02-10 (early Terra commissioning gap) - the animation effectively
#     starts 2000-02-18, 454 frames total. Instead of warping the grid through
#     terra::project() every frame, every map pixel is inverse-projected from
#     Robinson to lon/lat ONCE with PROJ (sf_project); since the projection never
#     changes, each frame is then just a lookup into that month's grid. Each frame
#     is written to its own PNG and already-rendered frames are skipped on rerun,
#     so the render can be interrupted and resumed; frames are stitched with av
#     once all are present, then the frame folder is deleted.

library(terra)
library(ncdf4)
library(sf)
library(rnaturalearth)
library(av)

# ---- Config -------------------------------------------------------------
DATA_PATH   <- "/Volumes/CLIMATEDATA/gpp_2000_2009.nc"
OUTPUT_PATH <- "/Volumes/CLIMATEDATA/gpp_2000_2009_robinson_R.mp4"
FRAMES_DIR  <- "/Volumes/CLIMATEDATA/_frames_tmp_gpp_robinson_R"
VAR_NAME    <- "gpp"

PROJ_CRS   <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"  # = ccrs.Robinson()
LONLAT_CRS <- "+proj=longlat +datum=WGS84 +no_defs"

FPS <- 4.7  # matches the original UW clip: ~454 frames / 98 seconds

# Test mode: render only the first MAX_FRAMES frames into a short "_test" MP4.
# Set to NULL for the full render (test frames are reused).
MAX_FRAMES <- NULL

VMIN  <- 0
VMAX  <- 0.12  # matches the observed value range from spot-checking
GAMMA <- 0.6   # <1 darkens mid-range values, pulling more of the colormap's
# saturated greens into the visible range instead of clustering
# everything in the pale end

# matplotlib "YlGn" (ColorBrewer 9-class, linearly interpolated to 256 colors).
N_COLORS <- 256
PALETTE <- colorRampPalette(c(
  "#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679",
  "#41ab5d", "#238443", "#006837", "#004529"
))(N_COLORS)

OCEAN_COLOR     <- "#0a1a3c"  # dark navy, matches the original clip's ocean
LAND_BASE_COLOR <- "#f2f2ea"  # off-white/gray base for non-vegetated land
COASTLINE_COLOR <- "#808080"  # matplotlib "gray" (R's "gray" is much lighter)

# Figure: 12 x 6 in at 150 dpi, same as the Python figure (even dims for H.264).
FIG_W_PX <- 1800
FIG_H_PX <- 900
FIG_DPI  <- 150
PT_Y <- (FIG_DPI / 72) / FIG_H_PX  # one typographic point, in figure fractions

# GPP colorbar placement/style. Figure-fraction [left, bottom, width, height].
COLORBAR_RECT       <- c(0.28, 0.085, 0.44, 0.035)
COLORBAR_LABEL      <- expression("GPP (kg C/m"^2 * "/8-day)")
COLORBAR_TICKS      <- c(0, 0.03, 0.06, 0.09, 0.12)
COLORBAR_TEXT_COLOR <- "black"

# Map axes placement, figure-fraction [left, bottom, width, height] -- same
# as the Python MAP_RECT. Like Cartopy, the map keeps Robinson's true aspect
# ratio, so it is height-limited and centered horizontally inside this box.
MAP_RECT <- c(0.01, 0.13, 0.98, 0.80)

# R's lwd = 1 is 1/96 in (0.75 pt); convert matplotlib point widths.
.pt_lwd <- function(pt) pt / 0.75


# ---- Data loading -------------------------------------------------------

.read_time <- function(nc) {
  vals  <- ncvar_get(nc, "time")
  units <- ncatt_get(nc, "time", "units")$value
  parts <- regmatches(units, regexec("^\\s*(\\w+)\\s+since\\s+(.+?)\\s*$", units))[[1]]
  if (length(parts) != 3) stop("Unrecognized time units: ", units)
  secs <- c(second = 1, seconds = 1, minute = 60, minutes = 60,
            hour = 3600, hours = 3600, day = 86400, days = 86400)[tolower(parts[2])]
  if (is.na(secs)) stop("Unsupported time step in units: ", units)
  origin <- as.POSIXct(sub("T", " ", sub("Z$", "", parts[3])), tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M", "%Y-%m-%d"))
  as.Date(origin + vals * secs)
}

# Reads one time step as an [x, y] matrix -- only a single frame is ever in
# memory (the R equivalent of chunks={"time": 1} + cache_frame_data=False).
.read_slice <- function(nc, t_idx) {
  dn <- vapply(nc$var[[VAR_NAME]]$dim, function(d) d$name, character(1))
  start <- rep(1, length(dn))
  count <- rep(-1, length(dn))
  start[dn == "time"] <- t_idx
  count[dn == "time"] <- 1
  m <- ncvar_get(nc, VAR_NAME, start = start, count = count, collapse_degen = TRUE)
  if (which(dn == "x") > which(dn == "y")) m <- t(m)
  m
}

# Project convention: force a real data read, not just a metadata open.
.verify_readable <- function(nc) {
  if (is.null(nc$var[[VAR_NAME]])) stop("Variable ", VAR_NAME, " missing from ", DATA_PATH)
  m <- .read_slice(nc, 1)
  if (all(is.na(m))) stop("First ", VAR_NAME, " slice in ", DATA_PATH, " is entirely NA")
  rm(m)
  invisible(gc())
}


# ---- Projection setup (done once) ------------------------------------------

.forward <- function(lon, lat) {
  sf::sf_project(LONLAT_CRS, PROJ_CRS, cbind(lon, lat), keep = TRUE, warn = FALSE)
}

# Map box in figure fractions and pixels, sized from PROJ's own Robinson extent.
.map_layout <- function() {
  x_max <- .forward(180, 0)[1, 1]
  y_max <- .forward(0, 90)[1, 2]
  px_h  <- round(MAP_RECT[4] * FIG_H_PX)
  px_w  <- round(px_h * x_max / y_max)
  w     <- px_w / FIG_W_PX
  left  <- MAP_RECT[1] + (MAP_RECT[3] - w) / 2
  list(x_max = x_max, y_max = y_max, px_w = px_w, px_h = px_h,
       left = left, right = left + w,
       bottom = MAP_RECT[2], top = MAP_RECT[2] + MAP_RECT[4])
}

# Projected meters -> figure fractions.
.to_fig <- function(X, Y, lay) {
  list(x = lay$left + (X + lay$x_max) / (2 * lay$x_max) * (lay$right - lay$left),
       y = lay$bottom + (Y + lay$y_max) / (2 * lay$y_max) * (lay$top - lay$bottom))
}

# Inverse-projects every map pixel center to lon/lat. Pixels outside the
# Robinson outline come back non-finite from PROJ and are left off the globe.
.pixel_lonlat <- function(lay) {
  xc <- ((seq_len(lay$px_w) - 0.5) / lay$px_w * 2 - 1) * lay$x_max
  yc <- (1 - (seq_len(lay$px_h) - 0.5) / lay$px_h * 2) * lay$y_max  # top row first
  X <- rep(xc, times = lay$px_h)                                   # row-major
  Y <- rep(yc, each = lay$px_w)
  ll <- sf::sf_project(PROJ_CRS, LONLAT_CRS, cbind(X, Y), keep = TRUE, warn = FALSE)
  inside <- is.finite(ll[, 1]) & is.finite(ll[, 2]) &
    abs(ll[, 1]) <= 180 & abs(ll[, 2]) <= 90
  list(lon = ll[inside, 1], lat = ll[inside, 2], inside = inside)
}

# Nearest data cell for each on-map pixel (NA where the pixel falls outside
# the data grid, e.g. beyond its latitude coverage -- never edge-smeared).
.data_index <- function(x, y, lon, lat) {
  if (max(x) > 180) lon <- lon %% 360
  nx <- length(x); ny <- length(y)
  dx <- x[2] - x[1]; dy <- y[2] - y[1]
  i <- round((lon - x[1]) / dx) + 1
  j <- round((lat - y[1]) / dy) + 1
  if (abs(nx * dx) > 359.9) i <- ((i - 1) %% nx) + 1   # global grid: wrap the seam
  ok <- i >= 1 & i <= nx & j >= 1 & j <= ny
  idx <- rep(NA_integer_, length(lon))
  idx[ok] <- as.integer(i[ok] + (j[ok] - 1) * nx)       # linear index into [x, y]
  idx
}

# Land/ocean base (Natural Earth 110m = Cartopy's default LAND feature),
# rasterized to 0.1 deg and sampled once at every on-map pixel.
.land_flags <- function(lon, lat) {
  land <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  template <- rast(nrows = 1800, ncols = 3600, xmin = -180, xmax = 180,
                   ymin = -90, ymax = 90, crs = "EPSG:4326")
  r <- rasterize(vect(land), template, field = 1, background = 0)
  m <- as.matrix(r, wide = TRUE)                     # [row = lat desc, col = lon]
  col <- pmin(pmax(floor((lon + 180) / 0.1) + 1, 1), 3600)
  row <- pmin(pmax(floor((90 - lat) / 0.1) + 1, 1), 1800)
  m[cbind(row, col)] == 1
}

# Coastlines (Natural Earth 110m) in figure fractions, NA-separated lines.
.coastline_fig <- function(lay) {
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
  xy  <- sf::st_coordinates(sf::st_transform(coast, PROJ_CRS))
  grp <- do.call(paste, as.data.frame(xy[, grepl("^L", colnames(xy)), drop = FALSE]))
  brk <- c(FALSE, grp[-1] != grp[-length(grp)])
  pos <- seq_len(nrow(xy)) + cumsum(brk)
  X <- Y <- rep(NA_real_, nrow(xy) + sum(brk))
  X[pos] <- xy[, 1]
  Y[pos] <- xy[, 2]
  .to_fig(X, Y, lay)
}

# Robinson outline (the map "spine"), traced along the +/-180 meridians and poles.
.outline_fig <- function(lay) {
  lat <- seq(-90, 90, length.out = 181)
  lon <- c(rep(-180, 181), seq(-180, 180, length.out = 181), rep(180, 181),
           seq(180, -180, length.out = 181))
  la  <- c(lat, rep(90, 181), rev(lat), rep(-90, 181))
  p <- .forward(lon, la)
  .to_fig(p[, 1], p[, 2], lay)
}


# ---- Rendering ----------------------------------------------------------

.format_date <- function(d) {
  # "Feb 18, 2000" -- month.abb is locale-independent, and as.integer drops the
  # leading zero portably (the Python %-d code is not portable).
  sprintf("%s %d, %s", month.abb[as.integer(format(d, "%m"))],
          as.integer(format(d, "%d")), format(d, "%Y"))
}

# PowerNorm + colormap: u = clip((v - VMIN) / (VMAX - VMIN), 0, 1)^GAMMA, split
# into N_COLORS equal bins. Out-of-range values take the end colors (matplotlib
# under/over default); NA returns NA so the base layer shows through.
.gpp_colors <- function(v) {
  out <- rep(NA_character_, length(v))
  ok  <- !is.na(v)
  u   <- pmin(pmax((v[ok] - VMIN) / (VMAX - VMIN), 0), 1)^GAMMA
  out[ok] <- PALETTE[pmin(floor(u * N_COLORS) + 1, N_COLORS)]
  out
}

.draw_gpp_colorbar <- function() {
  # Built directly from PALETTE / VMIN / VMAX / GAMMA, so it can't drift out of
  # sync with the map coloring. Like matplotlib's colorbar for a PowerNorm, the
  # colors are evenly spaced along the bar (true cell edges, N_COLORS + 1) and
  # the ticks sit at their power-scaled positions.
  l <- COLORBAR_RECT[1]; b <- COLORBAR_RECT[2]
  w <- COLORBAR_RECT[3]; h <- COLORBAR_RECT[4]
  xpos <- function(v) l + ((v - VMIN) / (VMAX - VMIN))^GAMMA * w
  
  edges <- seq(l, l + w, length.out = N_COLORS + 1)
  rect(edges[-(N_COLORS + 1)], b, edges[-1], b + h, col = PALETTE, border = NA)
  rect(l, b, l + w, b + h, border = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.6))
  
  segments(xpos(COLORBAR_TICKS), b, xpos(COLORBAR_TICKS), b - 3.5 * PT_Y,
           col = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.8))
  
  tick_label_top <- b - (3.5 + 3.5) * PT_Y   # tick length + matplotlib's default tick pad
  text(xpos(COLORBAR_TICKS), tick_label_top, labels = sprintf("%.2f", COLORBAR_TICKS),
       adj = c(0.5, 1), cex = 9 / 12, col = COLORBAR_TEXT_COLOR)
  text(l + w / 2, tick_label_top - (9 * 1.2 + 6) * PT_Y, labels = COLORBAR_LABEL,
       adj = c(0.5, 1), cex = 11 / 12, col = COLORBAR_TEXT_COLOR)
}

.render_frame <- function(frame_path, z, date_label, lay, geo) {
  # Base layer first (ocean/land), then GPP on top wherever it is not NA.
  cols <- geo$base_cols
  gpp  <- .gpp_colors(z[geo$data_idx])
  has  <- !is.na(gpp)
  cols[has] <- gpp[has]
  
  img <- rep("#FFFFFF", lay$px_w * lay$px_h)   # off-map area = figure background
  img[geo$inside] <- cols
  img <- as.raster(matrix(img, nrow = lay$px_h, byrow = TRUE))
  
  # Write to a partial file and rename on success, so an interrupted save can
  # never leave a truncated frame that the resume logic would then skip.
  tmp_path <- paste0(frame_path, ".partial")
  png(tmp_path, width = FIG_W_PX, height = FIG_H_PX, res = FIG_DPI, pointsize = 12, bg = "white")
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  rasterImage(img, lay$left, lay$bottom, lay$right, lay$top, interpolate = FALSE)
  lines(geo$coast$x, geo$coast$y, col = COASTLINE_COLOR, lwd = .pt_lwd(0.3))
  polygon(geo$outline$x, geo$outline$y, border = "black", lwd = .pt_lwd(0.8))
  
  title_y <- lay$top + 6 * PT_Y   # matplotlib's default title pad (6 pt)
  text(lay$left,  title_y, "Gross Primary Productivity", adj = c(0, 0), cex = 1)
  text(lay$right, title_y, date_label, adj = c(1, 0), cex = 1)
  
  .draw_gpp_colorbar()
  dev.off()
  file.rename(tmp_path, frame_path)
}

.list_frames <- function() {
  # ^ anchor also excludes AppleDouble "._" sidecars on the exFAT drive
  sort(list.files(FRAMES_DIR, pattern = "^frame_\\d{5}\\.png$", full.names = TRUE))
}

.assemble_video <- function(output_path, n_frames = NULL) {
  frame_files <- .list_frames()
  if (!is.null(n_frames)) frame_files <- head(frame_files, n_frames)
  message(sprintf("Assembling %d frames into %s at %s fps...", length(frame_files), output_path, FPS))
  if (file.exists(output_path)) file.remove(output_path)
  av::av_encode_video(frame_files, output = output_path, framerate = FPS,
                      codec = "libx264", verbose = FALSE)
  message("Saved animation to ", output_path)
}


# ---- Main ---------------------------------------------------------------

main <- function() {
  dir.create(FRAMES_DIR, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(FRAMES_DIR, pattern = "\\.partial$", full.names = TRUE))
  
  nc <- nc_open(DATA_PATH)
  on.exit(nc_close(nc), add = TRUE)
  .verify_readable(nc)
  
  dates    <- .read_time(nc)
  n_frames <- length(dates)
  
  message("Precomputing Robinson pixel lookup (one-time)...")
  lay <- .map_layout()
  pix <- .pixel_lonlat(lay)
  is_land <- .land_flags(pix$lon, pix$lat)
  geo <- list(
    inside    = pix$inside,
    data_idx  = .data_index(ncvar_get(nc, "x"), ncvar_get(nc, "y"), pix$lon, pix$lat),
    base_cols = ifelse(is_land, LAND_BASE_COLOR, OCEAN_COLOR),
    coast     = .coastline_fig(lay),
    outline   = .outline_fig(lay)
  )
  rm(pix, is_land)
  invisible(gc())
  
  test_mode     <- !is.null(MAX_FRAMES) && MAX_FRAMES < n_frames
  render_frames <- if (test_mode) MAX_FRAMES else n_frames
  if (test_mode) {
    message(sprintf("TEST MODE: rendering the first %d of %d frames.", render_frames, n_frames))
  } else {
    message(sprintf("Rendering %d frames at %s fps...", n_frames, FPS))
  }
  
  already_done <- length(.list_frames())
  if (already_done > 0) {
    message(sprintf("Resuming: %d frame(s) already rendered, skipping those.", already_done))
  }
  
  pb <- txtProgressBar(min = 0, max = render_frames, style = 3)
  for (i in seq_len(render_frames)) {
    setTxtProgressBar(pb, i)
    frame_path <- file.path(FRAMES_DIR, sprintf("frame_%05d.png", i - 1))
    if (file.exists(frame_path)) next
    
    z <- .read_slice(nc, i)
    .render_frame(frame_path, z, .format_date(dates[i]), lay, geo)
    rm(z)
    if (i %% 25 == 0) invisible(gc())
  }
  close(pb)
  
  if (test_mode) {
    .assemble_video(sub("\\.mp4$", "_test.mp4", OUTPUT_PATH), n_frames = render_frames)
    message("Set MAX_FRAMES <- NULL and rerun for the full render (test frames are reused).")
    return(invisible(NULL))
  }
  
  rendered_count <- length(.list_frames())
  if (rendered_count == n_frames) {
    .assemble_video(OUTPUT_PATH)
    unlink(FRAMES_DIR, recursive = TRUE)
  } else {
    message(sprintf(paste0("WARNING: expected %d frames but found %d -- skipping assembly. ",
                           "Rerun this script to fill in the missing frame(s) first."),
                    n_frames, rendered_count))
  }
}

main()