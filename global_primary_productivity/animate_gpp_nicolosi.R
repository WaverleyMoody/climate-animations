# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: animate_gpp_nicolosi.R
# Description: Generates the Gross Primary Productivity climatology
#     animation from MOD17A2HGF (2000-2009), rendered as a double-
#     hemisphere Nicolosi Globular projection.
# Note: For the Plate Carree, Robinson, and Foucaut projections, see the other
#     scripts in the gpp scripts folder. Data excludes 2000-01-01 through
#     2000-02-10 (early Terra commissioning gap) - the animation effectively
#     starts 2000-02-18, 454 frames total.
#
#     PROJ's `nicol` operation has no inverse transform, so terra::project()
#     (like Cartopy's GeoAxes in the Python version) can't be used. Instead,
#     following the lab's manual Nicolosi pipeline, a dense 0.1 deg lon/lat
#     grid is forward-projected ONCE per hemisphere and "rasterized" onto the
#     map pixels (each pixel records the lon/lat that landed on it); the few
#     pixels no sample hit are gap-filled with terra::focal(). That builds a
#     pixel -> lon/lat lookup table, so each frame is then just a lookup into
#     that month's grid. Natural Earth coastlines are clipped to each
#     hemisphere and forward-projected by hand with sf_project. Each frame is
#     written to its own PNG and already-rendered frames are skipped on rerun,
#     so the render can be interrupted and resumed; frames are stitched with av
#     once all are present, then the frame folder is deleted.

library(terra)
library(ncdf4)
library(sf)
library(rnaturalearth)
library(av)

sf::sf_use_s2(FALSE)   # planar lon/lat clipping to the hemisphere boxes

# ---- Config -------------------------------------------------------------
DATA_PATH   <- "/Volumes/CLIMATEDATA/gpp_2000_2009.nc"
OUTPUT_PATH <- "/Volumes/CLIMATEDATA/gpp_2000_2009_nicolosi_R.mp4"
FRAMES_DIR  <- "/Volumes/CLIMATEDATA/_frames_tmp_gpp_nicolosi_R"
VAR_NAME    <- "gpp"

WEST_LON <- -90.0   # hemisphere centers, matching the sibling Nicolosi scripts
EAST_LON <-  90.0
LONLAT_CRS <- "+proj=longlat +datum=WGS84 +no_defs"
.nicol_crs <- function(central_lon) sprintf("+proj=nicol +lon_0=%s +R=6371000", central_lon)

SAMPLE_STEP_DEG <- 0.1   # lon/lat sampling density for the forward-projected lookup

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
LINE_COLOR      <- "#808080"  # matplotlib "gray" (R's "gray" is much lighter)

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

# Hemisphere axes placement, figure-fraction [left, bottom, width, height] --
# same as the Python WEST_RECT/EAST_RECT. Each hemisphere disk keeps an equal
# aspect, so it is height-limited and centered horizontally inside its box.
WEST_RECT <- c(0.03, 0.13, 0.44, 0.80)
EAST_RECT <- c(0.53, 0.13, 0.44, 0.80)
TITLE_Y   <- 1.02  # title height as a fraction of the axes box, as in the Python version

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


# ---- Hemisphere setup (done once) ------------------------------------------

.forward <- function(lon, lat, central_lon) {
  sf::sf_project(LONLAT_CRS, .nicol_crs(central_lon), cbind(lon, lat), keep = TRUE, warn = FALSE)
}

# Disk placement in figure fractions and pixels.
.hemi_layout <- function(rect, central_lon) {
  edge_lats <- seq(-90, 90, length.out = 400)
  b1 <- .forward(rep(central_lon + 90, 400), edge_lats, central_lon)
  b2 <- .forward(rep(central_lon - 90, 400), edge_lats, central_lon)
  boundary <- rbind(b1, b2[400:1, ])                      # closed boundary circle
  r <- max(abs(boundary[, 1]), abs(boundary[, 2]), na.rm = TRUE)
  
  side_px <- round(rect[4] * FIG_H_PX)
  w <- side_px / FIG_W_PX
  left <- rect[1] + (rect[3] - w) / 2
  list(central_lon = central_lon, r = r, side_px = side_px, boundary = boundary,
       left = left, right = left + w, bottom = rect[2], top = rect[2] + rect[4],
       title_y = rect[2] + TITLE_Y * rect[4])
}

.to_fig <- function(X, Y, lay) {
  list(x = lay$left + (X + lay$r) / (2 * lay$r) * (lay$right - lay$left),
       y = lay$bottom + (Y + lay$r) / (2 * lay$r) * (lay$top - lay$bottom))
}

# Builds the pixel -> lon/lat lookup for one hemisphere by forward projection
# (nicol has no inverse). Returns lon/lat for every pixel inside the disk,
# row-major from the top-left, plus the inside mask.
.pixel_lonlat <- function(lay) {
  S  <- lay$side_px
  cl <- lay$central_lon
  
  offs <- seq(-90, 90, by = SAMPLE_STEP_DEG)             # lon offset from center
  lats <- seq(-90, 90, by = SAMPLE_STEP_DEG)
  off_s <- rep(offs, times = length(lats))
  lat_s <- rep(lats, each = length(offs))
  p <- .forward(cl + off_s, lat_s, cl)
  ok <- is.finite(p[, 1]) & is.finite(p[, 2])
  
  col <- pmin(pmax(floor((p[ok, 1] + lay$r) / (2 * lay$r) * S) + 1, 1), S)
  row <- pmin(pmax(floor((lay$r - p[ok, 2]) / (2 * lay$r) * S) + 1, 1), S)
  pidx <- (row - 1) * S + col                            # row-major pixel index
  
  off_px <- lat_px <- rep(NA_real_, S * S)
  off_px[pidx] <- off_s[ok]                              # a sample that landed here
  lat_px[pidx] <- lat_s[ok]
  rm(p, off_s, lat_s, ok, col, row, pidx)
  
  # Inside-the-disk mask from the pixel centers.
  cc <- ((seq_len(S) - 0.5) / S * 2 - 1)
  inside <- (rep(cc, times = S)^2 + rep(cc, each = S)^2) < 1
  
  # Gap-fill pixels inside the disk that no sample hit (focal mean of filled
  # neighbors, repeated until none remain) -- the lab's focal() gap-fill step.
  fill <- function(v) {
    r <- rast(matrix(v, nrow = S, byrow = TRUE))
    for (k in 1:20) {
      if (!any(is.na(values(r, mat = FALSE)[inside]))) break
      r <- focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    }
    as.vector(t(as.matrix(r, wide = TRUE)))
  }
  off_px <- fill(off_px)
  lat_px <- fill(lat_px)
  
  inside <- inside & !is.na(off_px) & !is.na(lat_px)
  lon <- ((cl + off_px[inside] + 180) %% 360) - 180      # back to -180..180
  list(lon = lon, lat = lat_px[inside], inside = inside)
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

# Land/ocean base (Natural Earth 110m, as in the Python version), rasterized
# to 0.1 deg once and sampled at every on-map pixel.
.build_land_mask <- function() {
  land <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  template <- rast(nrows = 1800, ncols = 3600, xmin = -180, xmax = 180,
                   ymin = -90, ymax = 90, crs = "EPSG:4326")
  as.matrix(rasterize(vect(land), template, field = 1, background = 0), wide = TRUE)
}
.land_flags <- function(mask, lon, lat) {
  col <- pmin(pmax(floor((lon + 180) / 0.1) + 1, 1), 3600)
  row <- pmin(pmax(floor((90 - lat) / 0.1) + 1, 1), 1800)
  mask[cbind(row, col)] == 1
}

# Coastlines clipped to the hemisphere box, then forward-projected by hand.
.coastline_fig <- function(coast, lay) {
  cl <- lay$central_lon
  bb <- sf::st_bbox(c(xmin = cl - 90, ymin = -90, xmax = cl + 90, ymax = 90),
                    crs = sf::st_crs(coast))
  geoms <- sf::st_geometry(suppressWarnings(sf::st_crop(coast, bb)))
  # Cropping can leave stray points where a line just touches the box edge.
  geoms <- geoms[sf::st_geometry_type(geoms) %in% c("LINESTRING", "MULTILINESTRING")]
  lines <- suppressWarnings(sf::st_cast(geoms, "LINESTRING"))
  xs <- ys <- numeric(0)
  for (ln in lines) {
    xy <- sf::st_coordinates(ln)
    p  <- .forward(xy[, 1], xy[, 2], cl)
    xs <- c(xs, p[, 1], NA)
    ys <- c(ys, p[, 2], NA)
  }
  .to_fig(xs, ys, lay)
}

.build_hemisphere <- function(rect, central_lon, x, y, land_mask, coast) {
  lay <- .hemi_layout(rect, central_lon)
  pix <- .pixel_lonlat(lay)
  is_land <- .land_flags(land_mask, pix$lon, pix$lat)
  list(
    lay       = lay,
    inside    = pix$inside,
    data_idx  = .data_index(x, y, pix$lon, pix$lat),
    base_cols = ifelse(is_land, LAND_BASE_COLOR, OCEAN_COLOR),
    coast     = .coastline_fig(coast, lay),
    boundary  = .to_fig(lay$boundary[, 1], lay$boundary[, 2], lay)
  )
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

.draw_hemisphere <- function(z, hemi) {
  # Base layer first (ocean/land), then GPP on top wherever it is not NA.
  cols <- hemi$base_cols
  gpp  <- .gpp_colors(z[hemi$data_idx])
  has  <- !is.na(gpp)
  cols[has] <- gpp[has]
  
  S <- hemi$lay$side_px
  img <- rep("#FFFFFF", S * S)   # outside the disk = figure background
  img[hemi$inside] <- cols
  img <- as.raster(matrix(img, nrow = S, byrow = TRUE))
  
  rasterImage(img, hemi$lay$left, hemi$lay$bottom, hemi$lay$right, hemi$lay$top,
              interpolate = FALSE)
  lines(hemi$coast$x, hemi$coast$y, col = LINE_COLOR, lwd = .pt_lwd(0.3))
  polygon(hemi$boundary$x, hemi$boundary$y, border = LINE_COLOR, lwd = .pt_lwd(0.6))
}

.render_frame <- function(frame_path, z, date_label, west, east) {
  # Write to a partial file and rename on success, so an interrupted save can
  # never leave a truncated frame that the resume logic would then skip.
  tmp_path <- paste0(frame_path, ".partial")
  png(tmp_path, width = FIG_W_PX, height = FIG_H_PX, res = FIG_DPI, pointsize = 12, bg = "white")
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  .draw_hemisphere(z, west)
  .draw_hemisphere(z, east)
  
  # Title y = TITLE_Y of the axes box, plus matplotlib's default 6 pt title pad.
  text(west$lay$left,  west$lay$title_y + 6 * PT_Y, "Gross Primary Productivity",
       adj = c(0, 0), cex = 1)
  text(east$lay$right, east$lay$title_y + 6 * PT_Y, date_label, adj = c(1, 0), cex = 1)
  
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
  x <- ncvar_get(nc, "x")
  y <- ncvar_get(nc, "y")
  
  land_mask <- .build_land_mask()
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
  
  message("Precomputing Western Hemisphere geometry...")
  west <- .build_hemisphere(WEST_RECT, WEST_LON, x, y, land_mask, coast)
  invisible(gc())
  message("Precomputing Eastern Hemisphere geometry...")
  east <- .build_hemisphere(EAST_RECT, EAST_LON, x, y, land_mask, coast)
  rm(land_mask)
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
    .render_frame(frame_path, z, .format_date(dates[i]), west, east)
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