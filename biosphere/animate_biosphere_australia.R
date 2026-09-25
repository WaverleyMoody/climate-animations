# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations
# Library by Professor John Michael Wallace.
#
# Script: animate_biosphere_australia.R
#
# Description: Generates the Global Biosphere animation from the merged land NDVI + ocean CHL dataset
#     (2000-2010, 0.05 deg grid). Two acts: (1) a full 360-degree world tour that plays through the
#     entire 2000-2010 monthly record while rotating, so the establishing shot shows real data
#     evolving; (2) a zoom into a fixed, close-up view of the focus region
#     (currently Australia) that replays the full 2000-2010 record a second time. Rendering is
#     resumable -- each frame is written to its own PNG and already-rendered frames are skipped on
#     rerun, so the render can be safely interrupted (e.g. closing the laptop lid) and picked up
#     again later without losing progress.
#
# Note: Land NDVI and ocean CHL are two independently-sourced variables in the same merged NetCDF,
#     each NaN over the other's domain, so they composite directly by layering (CHL first, NDVI on
#     top) with no additional masking logic needed. Residual CHL gaps (sea ice, quality-flag
#     exclusions) are real, expected features of ocean color data, not defects, and are rendered as
#     a neutral ocean base color rather than being filled in. The simulated satellite view uses the
#     PROJ near-sided perspective (nsper) math, the same projection behind Cartopy's
#     NearsidePerspective in the Python version. Instead of warping the 0.05 deg grid through
#     terra::project() every frame (26M cells, heavy on 8GB), each output pixel is inverse-projected
#     analytically to lon/lat and sampled nearest-neighbor from the month's grid; coastlines are
#     forward-projected and clipped at the horizon. For the Python version, see
#     animate_biosphere_australia.py.

library(terra)
library(ncdf4)
library(sf)
library(rnaturalearth)
library(av)

# --- Config -------------------------------------------------------------------

MERGED_DIR  <- "/Volumes/CLIMATEDATA/global_biosphere/merged"
FRAMES_DIR  <- "/Users/waverleymoody/CLIMATE ANIMATIONS/global_biosphere/frames_R"
OUTPUT_PATH <- "/Users/waverleymoody/CLIMATE ANIMATIONS/global_biosphere/animate_biosphere_australia_R.mp4"

# How many consecutive frames each month is held for before advancing to the
# next month. Camera motion (rotation, zoom) still updates every single frame
# for smoothness -- only the DATA changes at this slower rate.
HOLD_FRAMES_PER_MONTH <- 3

FPS <- 8   # each month is on screen for HOLD_FRAMES_PER_MONTH / FPS seconds

EASE_FRAMES <- 40   # transition frames between the world tour and the focus-region zoom

ROTATION_CENTRAL_LATITUDE <- 10.0   # slight tilt during the world tour

# Australia -- fixed, close-up focus region for the second act.
FOCUS_LON <- 135.0
FOCUS_LAT <- -25.0

# Simulated satellite altitude (meters). ~35.8 million m is geostationary
# altitude (whole-Earth disk, like GOES full-disk imagery) for the world tour;
# the focus-region zoom eases down to a much lower altitude.
WORLD_TOUR_SATELLITE_HEIGHT <- 35785831
FOCUS_SATELLITE_HEIGHT      <- 4000000

EARTH_RADIUS <- 6378137   # sphere radius (m) for the nsper math

# Figure: 10 x 10 in at 150 dpi, same as the Python figure (even dims for H.264).
FIG_PX  <- 1500
FIG_DPI <- 150
PT      <- (FIG_DPI / 72) / FIG_PX   # one typographic point, in figure-fraction units

# Map disk placement in figure fractions -- matches matplotlib's default
# single-subplot axes box ([0.125, 0.9] x [0.11, 0.88]), where the circular
# globe is limited by the axes height and centered in the box.
MAP_CX   <- 0.5125
MAP_CY   <- 0.495
MAP_SIDE <- 0.77
MAP_PX   <- round(MAP_SIDE * FIG_PX)   # pixel resolution of the globe image

N_COLORS <- 256   # matplotlib's default colormap resolution

# NDVI: matplotlib "YlGn" (ColorBrewer 9-class, linearly interpolated).
# NaN (ocean) is transparent so the ocean layer underneath shows through.
NDVI_PALETTE <- colorRampPalette(c(
  "#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679",
  "#41ab5d", "#238443", "#006837", "#004529"
))(N_COLORS)
NDVI_VMIN <- -0.2
NDVI_VMAX <- 1.0

# CHL: matplotlib "nipy_spectral" (its 21 segment stops at 0.05 spacing),
# log-scaled, since chlorophyll spans orders of magnitude (open-ocean deserts
# near 0.01 mg/m3 up to coastal blooms near 20+ mg/m3). NaN is transparent so
# the neutral ocean base color shows through.
NIPY_R <- c(0, .4667, .5333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .7333, .9333, 1, 1, 1, .8667, .8, .8)
NIPY_G <- c(0, 0, 0, 0, 0, .4667, .6, .6667, .6667, .6, .7333, .8667, 1, 1, .9333, .8, .6, 0, 0, 0, .8)
NIPY_B <- c(0, .5333, .6, .6667, .8667, .8667, .8667, .6667, .5333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .8)
CHL_PALETTE <- colorRampPalette(rgb(NIPY_R, NIPY_G, NIPY_B))(N_COLORS)
CHL_VMIN <- 0.01
CHL_VMAX <- 20.0

OCEAN_BASE_COLOR <- "#0a1a2f"   # neutral dark navy for ocean pixels with no CHL data
LAND_BASE_COLOR  <- "#3a3a3a"   # neutral gray for land pixels with no NDVI data (rare)
COASTLINE_COLOR  <- "#888888"

# CHL colorbar placement/style. Figure-fraction [left, bottom, width, height],
# tucked into the bottom-left so it clears the title and date label.
COLORBAR_RECT       <- c(0.06, 0.06, 0.30, 0.02)
COLORBAR_LABEL      <- expression("Chlorophyll-a (mg/m"^3 * ")")
COLORBAR_TICKS      <- c(0.01, 0.1, 1, 10)
COLORBAR_TEXT_COLOR <- "white"

# R's lwd = 1 is 1/96 in (0.75 pt); convert matplotlib point widths.
.pt_lwd <- function(pt) pt / 0.75


# --- Data loading ---------------------------------------------------------------

.list_merged_files <- function() {
  # ^ anchor also excludes AppleDouble "._" sidecars on the exFAT drive
  sort(list.files(MERGED_DIR, pattern = "^merged_biosphere_.*\\.nc$", full.names = TRUE))
}

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

# One row per month across all yearly files, sorted by date (mirrors
# open_mfdataset(combine="by_coords")).
.build_month_table <- function(ncs) {
  tbl <- do.call(rbind, lapply(seq_along(ncs), function(k) {
    d <- .read_time(ncs[[k]])
    data.frame(file_id = k, local_idx = seq_along(d), date = d)
  }))
  tbl <- tbl[order(tbl$date), ]
  rownames(tbl) <- NULL
  tbl
}

# Regular-grid description used for nearest-cell lookups.
.grid_spec <- function(lon, lat) {
  list(nlon = length(lon), lon0 = lon[1], dlon = lon[2] - lon[1],
       nlat = length(lat), lat0 = lat[1], dlat = lat[2] - lat[1])
}

# Reads one month of one variable as a [lon, lat] matrix -- only a single
# time slice is ever in memory (the R equivalent of chunks={"time": 1}).
.read_slice <- function(nc, var, t_idx) {
  dn <- vapply(nc$var[[var]]$dim, function(d) d$name, character(1))
  start <- rep(1, length(dn))
  count <- rep(-1, length(dn))
  start[dn == "time"] <- t_idx
  count[dn == "time"] <- 1
  m <- ncvar_get(nc, var, start = start, count = count, collapse_degen = TRUE)
  if (which(dn == "lon") > which(dn == "lat")) m <- t(m)
  m
}

# Project convention: force a real data read, not just a metadata open.
.verify_readable <- function(nc) {
  for (var in c("NDVI", "CHL")) {
    if (is.null(nc$var[[var]])) stop("Variable ", var, " missing from ", nc$filename)
    m <- .read_slice(nc, var, 1)
    if (all(is.na(m))) stop("First ", var, " slice in ", nc$filename, " is entirely NA")
    rm(m)
  }
  invisible(gc())
}

# Land/ocean base layer (Natural Earth 50m), rasterized once to 0.1 deg.
.build_land_mask <- function() {
  land <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
  template <- rast(nrows = 1800, ncols = 3600, xmin = -180, xmax = 180,
                   ymin = -90, ymax = 90, crs = "EPSG:4326")
  r <- rasterize(vect(land), template, field = 1, background = 0)
  list(values = t(as.matrix(r, wide = TRUE)),   # [lon, lat], lat descending
       grid = list(nlon = 3600, lon0 = -179.95, dlon = 0.1,
                   nlat = 1800, lat0 = 89.95, dlat = -0.1))
}

# Natural Earth 50m coastlines as lon/lat vectors with NA breaks between lines.
.build_coastline <- function() {
  coast <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
  xy  <- sf::st_coordinates(coast)
  grp <- do.call(paste, as.data.frame(xy[, grepl("^L", colnames(xy)), drop = FALSE]))
  brk <- c(FALSE, grp[-1] != grp[-length(grp)])
  pos <- seq_len(nrow(xy)) + cumsum(brk)
  lon <- lat <- rep(NA_real_, nrow(xy) + sum(brk))
  lon[pos] <- xy[, 1]
  lat[pos] <- xy[, 2]
  list(lon = lon, lat = lat)
}


# --- Helpers --------------------------------------------------------------------

.smoothstep <- function(t) {
  # Standard smoothstep easing (3t^2 - 2t^3), matching the SST animation's easing convention.
  t <- pmin(pmax(t, 0), 1)
  t * t * (3 - 2 * t)
}

.format_date <- function(d) {
  # month.abb is locale-independent
  paste(month.abb[as.integer(format(d, "%m"))], format(d, "%Y"))
}

.map_colors <- function(v, vmin, vmax, palette, log = FALSE) {
  # Same semantics as a matplotlib Normalize/LogNorm + colormap: out-of-range
  # values clip to the end colors, NaN (and <= 0 under log) stay transparent.
  out <- rep(NA_character_, length(v))
  ok  <- !is.na(v)
  if (log) {
    ok <- ok & v > 0
    t  <- (log10(v[ok]) - log10(vmin)) / (log10(vmax) - log10(vmin))
  } else {
    t <- (v[ok] - vmin) / (vmax - vmin)
  }
  idx <- pmin(pmax(floor(t * length(palette)) + 1, 1), length(palette))
  out[ok] <- palette[idx]
  out
}

.grid_index <- function(grid, lon, lat) {
  i <- (round(((lon - grid$lon0) %% 360) / grid$dlon) %% grid$nlon) + 1
  j <- pmin(pmax(round((lat - grid$lat0) / grid$dlat) + 1, 1), grid$nlat)
  cbind(i, j)
}

# Visible-disk radius (m) of the nsper projection: R * sqrt(h / (h + 2R)).
.disk_radius <- function(h) EARTH_RADIUS * sqrt(h / (h + 2 * EARTH_RADIUS))

# Near-sided perspective, spherical forward (Snyder 1987, eq. 23-1 to 23-5).
# Points beyond the horizon come back NA, which also breaks coastline lines there.
.nsper_forward <- function(lon, lat, lon0, lat0, h) {
  R <- EARTH_RADIUS
  P <- (R + h) / R
  phi  <- lat * pi / 180
  phi0 <- lat0 * pi / 180
  lam  <- (lon - lon0) * pi / 180
  cosc <- sin(phi0) * sin(phi) + cos(phi0) * cos(phi) * cos(lam)
  k <- (P - 1) / (P - cosc)
  x <- R * k * cos(phi) * sin(lam)
  y <- R * k * (cos(phi0) * sin(phi) - sin(phi0) * cos(phi) * cos(lam))
  hidden <- is.na(cosc) | cosc < 1 / P
  x[hidden] <- NA
  y[hidden] <- NA
  list(x = x, y = y)
}

# Near-sided perspective, spherical inverse (Snyder 1987, eq. 23-18 and 20-14/15).
.nsper_inverse <- function(x, y, lon0, lat0, h) {
  R <- EARTH_RADIUS
  P <- (R + h) / R
  phi0 <- lat0 * pi / 180
  rho  <- sqrt(x^2 + y^2)
  disc <- 1 - rho^2 * (P + 1) / (R^2 * (P - 1))
  disc[disc < 0] <- NA
  sinc <- (P - sqrt(disc)) / (R * (P - 1) / rho + rho / (R * (P - 1)))
  sinc[rho == 0] <- 0
  cosc <- sqrt(1 - sinc^2)
  rho_safe <- ifelse(rho == 0, 1, rho)
  phi <- asin(pmin(pmax(cosc * sin(phi0) + y * sinc * cos(phi0) / rho_safe, -1), 1))
  lam <- atan2(x * sinc, rho * cos(phi0) * cosc - y * sin(phi0) * sinc)
  list(lon = lon0 + lam * 180 / pi, lat = phi * 180 / pi)
}

# Normalized pixel-center coordinates of the square globe image, row-major
# from the top-left, restricted to the unit disk. Constant across frames.
.pixel_grid <- function(n) {
  centers <- ((seq_len(n) - 0.5) / n) * 2 - 1
  u <- rep(centers, times = n)
  v <- rep(rev(centers), each = n)
  inside <- u^2 + v^2 < 1
  list(u = u[inside], v = v[inside], inside = inside)
}

.frame_plan <- function(n_months) {
  # Three acts:
  #   1. World tour: one full 360-degree rotation while playing through every
  #      month in the record (n_months * HOLD_FRAMES_PER_MONTH frames).
  #   2. Ease: camera pans and zooms from wherever the tour ended to the fixed,
  #      close-up focus-region view (EASE_FRAMES long, data held on month 1 so
  #      the second act starts the record over from the beginning).
  #   3. Focus zoom: fixed close-up view, replaying the full record again.
  tour_n <- n_months * HOLD_FRAMES_PER_MONTH
  f <- 0:(tour_n - 1)

  tour <- data.frame(
    lon = 180 - 360 * (f / tour_n), lat = ROTATION_CENTRAL_LATITUDE,
    height = WORLD_TOUR_SATELLITE_HEIGHT,
    month_index = f %/% HOLD_FRAMES_PER_MONTH + 1, show_date = TRUE
  )

  tour_end_lon <- tour$lon[tour_n]
  t <- .smoothstep((0:(EASE_FRAMES - 1)) / EASE_FRAMES)
  ease <- data.frame(
    lon = tour_end_lon + t * (FOCUS_LON - tour_end_lon),
    lat = ROTATION_CENTRAL_LATITUDE + t * (FOCUS_LAT - ROTATION_CENTRAL_LATITUDE),
    height = WORLD_TOUR_SATELLITE_HEIGHT + t * (FOCUS_SATELLITE_HEIGHT - WORLD_TOUR_SATELLITE_HEIGHT),
    month_index = 1, show_date = FALSE
  )

  focus <- data.frame(
    lon = FOCUS_LON, lat = FOCUS_LAT, height = FOCUS_SATELLITE_HEIGHT,
    month_index = f %/% HOLD_FRAMES_PER_MONTH + 1, show_date = TRUE
  )

  rbind(tour, ease, focus)
}

.same_spec <- function(a, b) {
  a$lon == b$lon && a$lat == b$lat && a$height == b$height &&
    a$month_index == b$month_index && a$show_date == b$show_date
}

.list_frames <- function(frames_dir = FRAMES_DIR) {
  sort(list.files(frames_dir, pattern = "^frame_\\d{5}\\.png$", full.names = TRUE))
}


# --- Rendering ---------------------------------------------------------------

.draw_chl_colorbar <- function() {
  # Built directly from CHL_PALETTE / CHL_VMIN / CHL_VMAX, so it can't drift out
  # of sync with the ocean coloring. Cells are drawn with true edges
  # (N_COLORS + 1 breaks), not midpoints, so no color overflows the border.
  l <- COLORBAR_RECT[1]; b <- COLORBAR_RECT[2]
  w <- COLORBAR_RECT[3]; h <- COLORBAR_RECT[4]
  lo <- log10(CHL_VMIN); hi <- log10(CHL_VMAX)
  xpos <- function(v) l + (log10(v) - lo) / (hi - lo) * w

  edges <- seq(l, l + w, length.out = N_COLORS + 1)
  rect(edges[-(N_COLORS + 1)], b, edges[-1], b + h, col = CHL_PALETTE, border = NA)
  rect(l, b, l + w, b + h, border = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.6))

  # Minor ticks at 2-9 x each decade (LogLocator(subs="auto")), major at COLORBAR_TICKS.
  minor <- as.vector(outer(2:9, 10^(floor(lo):ceiling(hi))))
  minor <- minor[minor >= CHL_VMIN & minor <= CHL_VMAX]
  segments(xpos(minor), b, xpos(minor), b - 2.0 * PT,
           col = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.6))
  segments(xpos(COLORBAR_TICKS), b, xpos(COLORBAR_TICKS), b - 3.5 * PT,
           col = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.8))

  tick_label_top <- b - (3.5 + 3.5) * PT   # tick length + matplotlib's default tick pad
  text(xpos(COLORBAR_TICKS), tick_label_top, labels = format(COLORBAR_TICKS, drop0trailing = TRUE, trim = TRUE),
       adj = c(0.5, 1), cex = 8 / 12, col = COLORBAR_TEXT_COLOR)
  text(l + w / 2, tick_label_top - (8 * 1.2 + 4) * PT, labels = COLORBAR_LABEL,
       adj = c(0.5, 1), cex = 9 / 12, col = COLORBAR_TEXT_COLOR)
}

.render_frame <- function(frame_path, spec, ndvi, chl, data_grid, date_label, land, coast, pix) {
  rho_max <- .disk_radius(spec$height)

  # Inverse-project every on-disk pixel to lon/lat, then sample the grids.
  ll <- .nsper_inverse(pix$u * rho_max, pix$v * rho_max, spec$lon, spec$lat, spec$height)
  ok <- !is.na(ll$lat)
  disk_cols <- rep(OCEAN_BASE_COLOR, length(ll$lat))
  sub <- disk_cols[ok]

  # Layering mirrors the Python zorder: ocean/land base -> CHL -> NDVI.
  land_idx <- .grid_index(land$grid, ll$lon[ok], ll$lat[ok])
  sub[land$values[land_idx] == 1] <- LAND_BASE_COLOR

  data_idx <- .grid_index(data_grid, ll$lon[ok], ll$lat[ok])
  chl_cols <- .map_colors(chl[data_idx], CHL_VMIN, CHL_VMAX, CHL_PALETTE, log = TRUE)
  sub[!is.na(chl_cols)] <- chl_cols[!is.na(chl_cols)]
  ndvi_cols <- .map_colors(ndvi[data_idx], NDVI_VMIN, NDVI_VMAX, NDVI_PALETTE)
  sub[!is.na(ndvi_cols)] <- ndvi_cols[!is.na(ndvi_cols)]
  disk_cols[ok] <- sub

  img <- rep("#000000", MAP_PX^2)
  img[pix$inside] <- disk_cols
  img <- as.raster(matrix(img, nrow = MAP_PX, byrow = TRUE))

  # Write to a partial file and rename on success, so an interrupted save can
  # never leave a truncated frame that the resume logic would then skip.
  tmp_path <- paste0(frame_path, ".partial")
  png(tmp_path, width = FIG_PX, height = FIG_PX, res = FIG_DPI, pointsize = 12, bg = "black")
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), bg = "black", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  half <- MAP_SIDE / 2
  rasterImage(img, MAP_CX - half, MAP_CY - half, MAP_CX + half, MAP_CY + half,
              interpolate = FALSE)

  cl <- .nsper_forward(coast$lon, coast$lat, spec$lon, spec$lat, spec$height)
  lines(MAP_CX + cl$x / rho_max * half, MAP_CY + cl$y / rho_max * half,
        col = COASTLINE_COLOR, lwd = .pt_lwd(0.4))

  text(0.02, 0.96, "Global Biosphere", adj = c(0, 1), cex = 16 / 12, font = 2, col = "white")
  if (!is.null(date_label)) {
    text(0.98, 0.96, date_label, adj = c(1, 1), cex = 14 / 12, col = "white")
  }

  .draw_chl_colorbar()
  dev.off()
  file.rename(tmp_path, frame_path)
}

.assemble_video <- function(frames_dir, output_path, fps) {
  # Stitches the rendered PNG frames into the final MP4 (libx264, yuv420p via av).
  frame_files <- .list_frames(frames_dir)
  if (length(frame_files) == 0) {
    message("No frames found in ", frames_dir, " -- nothing to assemble.")
    return(invisible(NULL))
  }
  message(sprintf("Assembling %d frames into %s at %d fps...", length(frame_files), output_path, fps))
  if (file.exists(output_path)) file.remove(output_path)
  av::av_encode_video(frame_files, output = output_path, framerate = fps,
                      codec = "libx264", verbose = FALSE)
  message("Saved video to ", output_path)
}


# --- Main ---------------------------------------------------------------------

main <- function() {
  dir.create(FRAMES_DIR, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(FRAMES_DIR, pattern = "\\.partial$", full.names = TRUE))

  files <- .list_merged_files()
  if (length(files) == 0) stop("No merged_biosphere_*.nc files found in ", MERGED_DIR)
  ncs <- lapply(files, nc_open)
  on.exit(invisible(lapply(ncs, nc_close)), add = TRUE)

  .verify_readable(ncs[[1]])
  months    <- .build_month_table(ncs)
  n_months  <- nrow(months)
  data_grid <- .grid_spec(ncvar_get(ncs[[1]], "lon"), ncvar_get(ncs[[1]], "lat"))

  plan <- .frame_plan(n_months)
  total_frames <- nrow(plan)

  land  <- .build_land_mask()
  coast <- .build_coastline()
  pix   <- .pixel_grid(MAP_PX)

  already_done <- length(.list_frames())
  if (already_done > 0) {
    message(sprintf("Resuming: %d frame(s) already rendered, skipping those.", already_done))
  }

  message("Rendering Global Biosphere frames")
  pb <- txtProgressBar(min = 0, max = total_frames, style = 3)
  cached_month <- -1
  ndvi <- chl <- NULL

  for (k in seq_len(total_frames)) {
    setTxtProgressBar(pb, k)
    frame_path <- file.path(FRAMES_DIR, sprintf("frame_%05d.png", k - 1))
    if (file.exists(frame_path)) next   # already rendered in a previous run

    spec <- plan[k, ]

    # Focus-act frames that hold the same month with a fixed camera are
    # pixel-identical to the previous frame -- copy instead of re-rendering.
    if (k > 1 && .same_spec(spec, plan[k - 1, ])) {
      prev_path <- file.path(FRAMES_DIR, sprintf("frame_%05d.png", k - 2))
      if (file.exists(prev_path)) {
        tmp_path <- paste0(frame_path, ".partial")
        file.copy(prev_path, tmp_path, overwrite = TRUE)
        file.rename(tmp_path, frame_path)
        next
      }
    }

    m <- spec$month_index
    if (m != cached_month) {
      ndvi <- chl <- NULL
      invisible(gc())
      nc   <- ncs[[months$file_id[m]]]
      ndvi <- .read_slice(nc, "NDVI", months$local_idx[m])
      chl  <- .read_slice(nc, "CHL",  months$local_idx[m])
      cached_month <- m
    }

    date_label <- if (spec$show_date) .format_date(months$date[m]) else NULL
    .render_frame(frame_path, spec, ndvi, chl, data_grid, date_label, land, coast, pix)
  }
  close(pb)

  rendered_count <- length(.list_frames())
  message(sprintf("All %d frames rendered to %s", total_frames, FRAMES_DIR))
  message(sprintf("Estimated runtime at %d fps: %.0f seconds", FPS, total_frames / FPS))

  if (rendered_count == total_frames) {
    .assemble_video(FRAMES_DIR, OUTPUT_PATH, FPS)
  } else {
    message(sprintf(paste0("WARNING: expected %d frames but found %d -- skipping assembly. ",
                           "Rerun this script to fill in the missing frame(s) first."),
                    total_frames, rendered_count))
  }
}

if (sys.nframe() == 0) main()