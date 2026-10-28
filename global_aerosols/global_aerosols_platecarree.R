# =============================================================================
# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations
# Library by Professor John Michael Wallace.
#
# Script: global_aerosols_platecarree.R
# Description: Generates the global aerosols animation from NASA MERRA-2 hourly
#     aerosol extinction AOT (dust, carbonaceous, sulfate, sea salt) and total
#     precipitable water vapor (2018-08-01 to 2020-01-31), rendered in the
#     PlateCarree projection.
# Note: DATA_DIR can hold either the daily hourly files from
#     download_aerosols_merra2.py (merra2_aerosols_YYYYMMDD.nc) or the monthly
#     3-hourly release files (global_aerosols_merra2_3hourly_YYYY-MM.nc); the
#     rendered frames are identical.
#
# Run:
#     # 1) Set PREVIEW_ONLY <- TRUE, check the PNG
#     # 2) Set PREVIEW_ONLY <- FALSE for the full render
#     caffeinate -i Rscript global_aerosols_platecarree.R
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(av)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
DATA_DIR        <- "/Volumes/CLIMATEDATA/global_aerosols/merra2_hourly"
BLUE_MARBLE_DIR <- "/Volumes/CLIMATEDATA/blue_marble/2004"

LOCAL_DIR <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "global_aerosols", "_render")  # local APFS
FINAL_DIR <- "/Volumes/CLIMATEDATA/global_aerosols/animations"                              # exFAT
NE_DIR    <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "global_aerosols", "_naturalearth")
OUT_NAME  <- "global_aerosols_2018_2020_platecarree_R.mp4"

PREVIEW_ONLY <- FALSE                      # TRUE: save one PNG and stop
PREVIEW_TIME <- "2019-08-02 14:30"        # matches the reference frame (UTC)

HOUR_STEP <- 3                            # render 00:30, 03:30, ... UTC -> 4,392 frames
FPS       <- 30                           # 4,392 frames -> ~2:26 (original 2:13)

# Composited map raster (2:1 equirectangular) and full video frame, in pixels
WIDTH  <- 1920; HEIGHT <- 960
FIG_W  <- 1920; FIG_H  <- 1200
DPI    <- 100

TITLE    <- "Global Aerosols, August 2018 – January 2020"
SUBTITLE <- "NASA MERRA-2 hourly aerosol optical thickness (550 nm) by species, with total precipitable water vapor"
CREDIT   <- "Data: NASA GMAO MERRA-2 (M2T1NXAER, M2T1NXSLV)  ·  Base map: NASA Blue Marble Next Generation (2004)"

LON_TICKS <- seq(-180, 180, by = 90)
LAT_TICKS <- seq(-90, 90, by = 45)
GRIDLINES <- TRUE

OCEAN_COLOR   <- "#0b2550"                # deep-ocean base under the sea salt layer
LAND_KEY_GRAY <- 0.55                     # backdrop used to draw the land-layer color keys
SEA_ICE_LUM   <- 0.35                     # Blue Marble ocean pixels brighter than this = sea ice

# Layers, composited bottom -> top.
#   vars:   summed MERRA-2 variables
#   scale:  "log" or "linear" between vmin and vmax
#   colors: color ramp from low to high
#   alpha:  maximum opacity; gamma shapes how fast opacity rises
#   ocean_only: draw only over ocean pixels
#   label, ticks, key_bg: color-key text, tick values, and backdrop ("ocean" or "land")
LAYERS <- list(
  list(name = "sea salt", vars = "SSEXTTAU", scale = "log", vmin = 0.02, vmax = 0.40,
       colors = c("#0d3a6e", "#2a7fb8", "#8fd4d9", "#f0eebf"),
       alpha = 0.95, gamma = 0.7, ocean_only = TRUE,
       label = "Sea salt AOT", ticks = c(0.02, 0.05, 0.1, 0.2, 0.4), key_bg = "ocean"),
  list(name = "water vapor", vars = "TQV", scale = "linear", vmin = 20, vmax = 65,
       colors = c("#c8ccd0", "#e6e8ea"),
       alpha = 0.25, gamma = 1.5, ocean_only = FALSE,
       label = expression("Water vapor (kg m"^-2 * ")"), ticks = c(20, 35, 50, 65), key_bg = "ocean"),
  list(name = "sulfate", vars = "SUEXTTAU", scale = "log", vmin = 0.05, vmax = 0.60,
       colors = c("#d9d9d9", "#ffffff"),
       alpha = 0.75, gamma = 1.0, ocean_only = FALSE,
       label = "Sulfate AOT", ticks = c(0.05, 0.1, 0.2, 0.6), key_bg = "land"),
  list(name = "carbon / smoke", vars = c("BCEXTTAU", "OCEXTTAU"), scale = "log", vmin = 0.04, vmax = 3.0,
       colors = c("#3f9e4a", "#7cc25a", "#c8331b", "#ff5a1a", "#ffd27a", "#fff6e0"),
       alpha = 0.95, gamma = 0.6, ocean_only = FALSE,
       label = "Smoke / carbon AOT", ticks = c(0.04, 0.1, 0.3, 1, 3), key_bg = "land"),
  list(name = "dust", vars = "DUEXTTAU", scale = "log", vmin = 0.05, vmax = 2.0,
       colors = c("#c56a1c", "#ef9a2e", "#f9c35a", "#fde9b0"),
       alpha = 0.92, gamma = 0.7, ocean_only = FALSE,
       label = "Dust AOT", ticks = c(0.05, 0.1, 0.2, 0.5, 1, 2), key_bg = "land")
)
KEY_ORDER <- c("dust", "carbon / smoke", "sulfate", "sea salt", "water vapor")   # left -> right
ALL_VARS  <- sort(unique(unlist(lapply(LAYERS, `[[`, "vars"))))

# Figure layout (fractions of the frame); map box is exactly 2:1 in pixels
MAP_X0 <- 0.07
MAP_W  <- 0.86
MAP_H  <- MAP_W * (FIG_W / FIG_H) / 2
MAP_Y0 <- 0.215
MAP_X1 <- MAP_X0 + MAP_W
MAP_Y1 <- MAP_Y0 + MAP_H

OCEAN_RGB <- as.numeric(col2rgb(OCEAN_COLOR)) / 255
PT        <- 12                           # device pointsize; cex = matplotlib fontsize / PT

# -----------------------------------------------------------------------------
# Target grid (pixel centers). Pixel vectors are column-major H x W:
# row 1 = north edge, column 1 = 180W.
# -----------------------------------------------------------------------------
LAT_C <- 90 - (seq_len(HEIGHT) - 0.5) * 180 / HEIGHT
LON_C <- -180 + (seq_len(WIDTH) - 0.5) * 360 / WIDTH

template_rast <- function() {
  terra::rast(nrows = HEIGHT, ncols = WIDTH, xmin = -180, xmax = 180, ymin = -90, ymax = 90,
              crs = "EPSG:4326")
}

terra_to_pixels <- function(v) {
  # terra cell order is row-major from the top-left; convert to column-major H x W
  as.vector(matrix(v, nrow = HEIGHT, ncol = WIDTH, byrow = TRUE))
}

# -----------------------------------------------------------------------------
# Static fields (built once)
# -----------------------------------------------------------------------------
build_interp <- function(src_lat, src_lon) {
  # Bilinear weights from the MERRA-2 grid to every output pixel; lon wraps at 180
  nlat <- length(src_lat); nlon <- length(src_lon)
  fr <- (LAT_C - src_lat[1]) / (src_lat[2] - src_lat[1])
  fr <- pmin(pmax(fr, 0), nlat - 1)
  r0 <- pmin(floor(fr), nlat - 2); wr <- fr - r0
  fc <- (LON_C - src_lon[1]) / (src_lon[2] - src_lon[1])
  c0 <- floor(fc); wc <- fc - c0
  c1 <- (c0 + 1) %% nlon; c0 <- c0 %% nlon
  
  R0 <- rep(r0, times = WIDTH); WR <- rep(wr, times = WIDTH)
  C0 <- rep(c0, each = HEIGHT); C1 <- rep(c1, each = HEIGHT); WC <- rep(wc, each = HEIGHT)
  
  # ncvar_get returns [lon, lat]: 1-based linear index = lon + lat * nlon + 1 (0-based lon/lat)
  list(i00 = C0 + R0 * nlon + 1,       i01 = C1 + R0 * nlon + 1,
       i10 = C0 + (R0 + 1) * nlon + 1, i11 = C1 + (R0 + 1) * nlon + 1,
       w00 = (1 - WR) * (1 - WC), w01 = (1 - WR) * WC,
       w10 = WR * (1 - WC),       w11 = WR * WC)
}

upsample <- function(f, ip) {
  f[ip$i00] * ip$w00 + f[ip$i01] * ip$w01 + f[ip$i10] * ip$w10 + f[ip$i11] * ip$w11
}

land_mask <- function() {
  cache <- file.path(LOCAL_DIR, sprintf("land_mask_%dx%d_R.rds", WIDTH, HEIGHT))
  if (file.exists(cache)) return(readRDS(cache))
  message("Building land mask (one-time)...")
  dir.create(NE_DIR, recursive = TRUE, showWarnings = FALSE)
  # Newer rnaturalearth saves a .gpkg, older versions a .shp — accept either
  find_land <- function() {
    hits <- file.path(NE_DIR, c("ne_50m_land.gpkg", "ne_50m_land.shp"))
    hits[file.exists(hits)][1]
  }
  if (is.na(find_land())) {
    ne_download(scale = 50, type = "land", category = "physical",
                destdir = NE_DIR, load = FALSE)
  }
  if (is.na(find_land())) stop("Natural Earth land file not found in ", NE_DIR)
  land <- st_read(find_land(), quiet = TRUE)
  r <- terra::rasterize(terra::vect(land), template_rast(), field = 1, background = 0)
  mask <- terra_to_pixels(terra::values(r)[, 1]) == 1
  saveRDS(mask, cache)
  mask
}

base_map <- function(month, land) {
  # Desaturated Blue Marble on land/sea ice, navy ocean
  f <- file.path(BLUE_MARBLE_DIR, sprintf("world.2004%02d.3x5400x2700.jpg", month))
  src <- suppressWarnings(terra::rast(f))
  # Rebuild in memory on a global lon/lat grid: resampling the un-georeferenced JPEG
  # directly from disk returns scattered NA cells
  r <- terra::rast(nrows = terra::nrow(src), ncols = terra::ncol(src), nlyrs = 3,
                   xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = "EPSG:4326")
  r <- terra::setValues(r, terra::values(src)[, 1:3])
  rm(src)
  rr <- terra::resample(r, template_rast(), method = "average")
  rm(r); invisible(gc())
  v <- terra::values(rr) / 255
  R <- terra_to_pixels(v[, 1]); G <- terra_to_pixels(v[, 2]); B <- terra_to_pixels(v[, 3])
  lum  <- 0.299 * R + 0.587 * G + 0.114 * B
  gray <- 0.35 + 0.65 * lum
  water <- (!land) & (lum < SEA_ICE_LUM)
  list(r = ifelse(water, OCEAN_RGB[1], gray),
       g = ifelse(water, OCEAN_RGB[2], gray),
       b = ifelse(water, OCEAN_RGB[3], gray))
}

ramp_lut <- function(colors) colorRamp(colors)(seq(0, 1, length.out = 256)) / 255

# -----------------------------------------------------------------------------
# Per-frame compositing
# -----------------------------------------------------------------------------
intensity <- function(x, L) {
  x <- pmax(x, 1e-6)
  t <- if (L$scale == "log") {
    (log10(x) - log10(L$vmin)) / (log10(L$vmax) - log10(L$vmin))
  } else {
    (x - L$vmin) / (L$vmax - L$vmin)
  }
  pmin(pmax(t, 0), 1)
}

composite <- function(fields, base, ocean, ip, luts) {
  r <- base$r; g <- base$g; b <- base$b
  for (k in seq_along(LAYERS)) {
    L <- LAYERS[[k]]; lut <- luts[[k]]
    raw <- Reduce(`+`, fields[L$vars])
    raw[is.na(raw)] <- 0
    t  <- intensity(upsample(raw, ip), L)
    ix <- as.integer(t * 255) + 1L
    a  <- L$alpha * t^L$gamma
    if (isTRUE(L$ocean_only)) a <- a * ocean
    r <- r * (1 - a) + lut[ix, 1] * a
    g <- g * (1 - a) + lut[ix, 2] * a
    b <- b * (1 - a) + lut[ix, 3] * a
  }
  array(pmin(pmax(c(r, g, b), 0), 1), dim = c(HEIGHT, WIDTH, 3))
}

# -----------------------------------------------------------------------------
# Drawing: map, axes labels, title, color keys
# -----------------------------------------------------------------------------
key_colors <- function(L, lut) {
  t  <- seq(0, 1, length.out = 256)
  a  <- L$alpha * t^L$gamma
  bg <- if (L$key_bg == "ocean") OCEAN_RGB else rep(LAND_KEY_GRAY, 3)
  rgb(bg[1] * (1 - a) + lut[, 1] * a,
      bg[2] * (1 - a) + lut[, 2] * a,
      bg[3] * (1 - a) + lut[, 3] * a)
}

lon_label <- function(v) {
  if (abs(v) == 180) return("180°")
  if (v == 0) return("0°")
  sprintf("%d°%s", as.integer(abs(v)), if (v > 0) "E" else "W")
}

lat_label <- function(v) {
  if (v == 0) return("0°")
  sprintf("%d°%s", as.integer(abs(v)), if (v > 0) "N" else "S")
}

lon_x <- function(lon) MAP_X0 + (lon + 180) / 360 * MAP_W
lat_y <- function(lat) MAP_Y0 + (lat + 90) / 180 * MAP_H

draw_frame <- function(img, stamp, keys) {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  # Map
  rasterImage(img, MAP_X0, MAP_Y0, MAP_X1, MAP_Y1, interpolate = TRUE)
  if (GRIDLINES) {
    grid_col <- rgb(1, 1, 1, 0.22)
    segments(lon_x(LON_TICKS), MAP_Y0, lon_x(LON_TICKS), MAP_Y1, col = grid_col, lwd = 0.8)
    segments(MAP_X0, lat_y(LAT_TICKS), MAP_X1, lat_y(LAT_TICKS), col = grid_col, lwd = 0.8)
  }
  rect(MAP_X0, MAP_Y0, MAP_X1, MAP_Y1, lwd = 1.1)
  
  # Ticks on all four sides, labels on the left and bottom only
  tx <- 5.5 / FIG_W; ty <- 5.5 / FIG_H
  segments(lon_x(LON_TICKS), MAP_Y0, lon_x(LON_TICKS), MAP_Y0 - ty, lwd = 1.1)
  segments(lon_x(LON_TICKS), MAP_Y1, lon_x(LON_TICKS), MAP_Y1 + ty, lwd = 1.1)
  segments(MAP_X0, lat_y(LAT_TICKS), MAP_X0 - tx, lat_y(LAT_TICKS), lwd = 1.1)
  segments(MAP_X1, lat_y(LAT_TICKS), MAP_X1 + tx, lat_y(LAT_TICKS), lwd = 1.1)
  text(lon_x(LON_TICKS), MAP_Y0 - ty - 0.006, sapply(LON_TICKS, lon_label),
       adj = c(0.5, 1), cex = 16 / PT)
  text(MAP_X0 - tx - 0.004, lat_y(LAT_TICKS), sapply(LAT_TICKS, lat_label),
       adj = c(1, 0.5), cex = 16 / PT)
  
  # Title, subtitle, timestamp
  text(MAP_X0, MAP_Y1 + 0.050, TITLE, adj = c(0, 0), cex = 32 / PT, font = 2)
  text(MAP_X0, MAP_Y1 + 0.018, SUBTITLE, adj = c(0, 0), cex = 18 / PT, col = "#333333")
  text(MAP_X1, MAP_Y1 + 0.050, stamp, adj = c(1, 0), cex = 28 / PT, family = "mono")
  
  # Color keys: one horizontal bar per layer
  n <- length(KEY_ORDER); gap <- 0.05
  bar_w <- (MAP_W - gap * (n - 1)) / n
  bar_y <- 0.105; bar_h <- 0.022
  for (i in seq_len(n)) {
    kk <- keys[[KEY_ORDER[i]]]; L <- kk$layer
    x0 <- MAP_X0 + (i - 1) * (bar_w + gap); x1 <- x0 + bar_w
    rasterImage(as.raster(matrix(kk$colors, nrow = 1)), x0, bar_y, x1, bar_y + bar_h,
                interpolate = FALSE)
    rect(x0, bar_y, x1, bar_y + bar_h, lwd = 0.8)
    pos <- if (L$scale == "log") {
      (log10(L$ticks) - log10(L$vmin)) / (log10(L$vmax) - log10(L$vmin))
    } else {
      (L$ticks - L$vmin) / (L$vmax - L$vmin)
    }
    tick_x <- x0 + pos * bar_w
    segments(tick_x, bar_y, tick_x, bar_y - ty, lwd = 0.8)
    text(tick_x, bar_y - ty - 0.005, as.character(L$ticks), adj = c(0.5, 1), cex = 15 / PT)
    label_y <- bar_y + bar_h + 0.009 - if (is.expression(L$label)) 0.004 else 0   # superscript raises plotmath
    text((x0 + x1) / 2, label_y, L$label, adj = c(0.5, 0), cex = 16 / PT)
  }
  
  text(MAP_X1, 0.015, CREDIT, adj = c(1, 0), cex = 13 / PT, col = "#555555")
}

# -----------------------------------------------------------------------------
# Files, time, verification
# -----------------------------------------------------------------------------
DAILY_PATTERN   <- "^merra2_aerosols_[0-9]{8}\\.nc$"                       # download script output
MONTHLY_PATTERN <- "^global_aerosols_merra2_3hourly_[0-9]{4}-[0-9]{2}\\.nc$"  # GitHub Release files

data_files <- function() {
  # Daily hourly files if present, otherwise monthly 3-hourly release files
  pick <- function(pattern) {
    f <- list.files(DATA_DIR, pattern = pattern, full.names = TRUE)
    sort(f[!startsWith(basename(f), "._")])
  }
  daily <- pick(DAILY_PATTERN)
  if (length(daily)) return(list(files = daily, kind = "daily"))
  monthly <- pick(MONTHLY_PATTERN)
  if (length(monthly)) return(list(files = monthly, kind = "monthly"))
  stop("No merra2_aerosols_YYYYMMDD.nc or global_aerosols_merra2_3hourly_YYYY-MM.nc files in ", DATA_DIR)
}

file_for <- function(t, kind) {
  if (kind == "daily") {
    file.path(DATA_DIR, sprintf("merra2_aerosols_%s.nc", format(t, "%Y%m%d", tz = "UTC")))
  } else {
    file.path(DATA_DIR, sprintf("global_aerosols_merra2_3hourly_%s.nc", format(t, "%Y-%m", tz = "UTC")))
  }
}

render_this <- function(t) as.integer(format(t, "%H", tz = "UTC")) %% HOUR_STEP == 0

read_times <- function(nc) {
  u      <- ncatt_get(nc, "time", "units")$value          # e.g. "minutes since 2019-08-02 00:30:00"
  unit   <- sub(" since.*", "", u)
  origin <- as.POSIXct(gsub("T", " ", sub(".*since ", "", u)), tz = "UTC")
  mult   <- c(seconds = 1, minutes = 60, hours = 3600, days = 86400)[[unit]]
  origin + as.numeric(ncvar_get(nc, "time")) * mult
}

open_data <- function(path) {
  # Daily files (24 steps) are read whole; monthly files frame by frame to limit memory.
  # get(i) returns a named list of [lon, lat] matrices for timestep i.
  nc    <- nc_open(path)
  times <- read_times(nc)
  out   <- list(times = times,
                lat   = as.numeric(ncvar_get(nc, "lat")),
                lon   = as.numeric(ncvar_get(nc, "lon")))
  if (length(times) <= 24) {
    arrs <- setNames(lapply(ALL_VARS, function(v) ncvar_get(nc, v)), ALL_VARS)  # [lon, lat, time]
    nc_close(nc)
    out$get   <- function(i) lapply(arrs, function(a) a[, , i])
    out$close <- function() invisible()
  } else {
    out$get <- function(i) {
      setNames(lapply(ALL_VARS, function(v) ncvar_get(nc, v, start = c(1, 1, i), count = c(-1, -1, 1))),
               ALL_VARS)
    }
    out$close <- function() nc_close(nc)
  }
  out
}

stamp_text <- function(t) format(t, "%Y-%m-%d %H:%M UTC", tz = "UTC")

verify_mp4 <- function(path, expected_frames) {
  if (!file.exists(path) || file.size(path) == 0) stop(basename(path), " is missing or empty")
  info <- av_media_info(path)
  # av's frame-count field can be off by one for H.264; duration x fps is exact
  n_frames <- as.integer(round(info$duration * FPS))
  if (n_frames != expected_frames) {
    stop(sprintf("%s: %d frames, expected %d", basename(path), n_frames, expected_frames))
  }
  message(sprintf("Verified %s: %d frames, %.1f s", basename(path), n_frames, info$duration))
}

move_to_final <- function(local_path) {
  # file.rename() fails across filesystems (APFS -> exFAT): copy, check size, delete
  dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
  final_path <- file.path(FINAL_DIR, basename(local_path))
  if (file.exists(final_path)) file.remove(final_path)
  ok <- file.copy(local_path, final_path)
  if (!ok || file.size(final_path) != file.size(local_path)) stop("Copy to ", FINAL_DIR, " failed")
  file.remove(local_path)
  final_path
}

open_device <- function(path, type = c("png", "jpeg")) {
  type <- match.arg(type)
  if (type == "png") {
    png(path, width = FIG_W, height = FIG_H, res = DPI, pointsize = PT, bg = "white")
  } else {
    jpeg(path, width = FIG_W, height = FIG_H, res = DPI, pointsize = PT, bg = "white", quality = 95)
  }
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  if (!dir.exists(DATA_DIR)) stop(DATA_DIR, " not found. Is /Volumes/CLIMATEDATA mounted?")
  dir.create(LOCAL_DIR, recursive = TRUE, showWarnings = FALSE)
  
  found <- data_files()
  files <- found$files; kind <- found$kind
  message(sprintf("Found %d %s files", length(files), kind))
  
  first <- open_data(files[1])
  first$close()
  stopifnot(isTRUE(all.equal(first$lon[1], -180)), isTRUE(all.equal(first$lat[1], -90)))
  ip <- build_interp(first$lat, first$lon)
  rm(first); invisible(gc())
  
  land  <- land_mask()
  ocean <- as.numeric(!land)
  luts  <- lapply(LAYERS, function(L) ramp_lut(L$colors))
  keys  <- setNames(lapply(seq_along(LAYERS), function(k) {
    list(layer = LAYERS[[k]], colors = key_colors(LAYERS[[k]], luts[[k]]))
  }), vapply(LAYERS, `[[`, "", "name"))
  
  # Only the current month's base map is kept in memory
  base_cache <- list(month = NA_integer_, base = NULL)
  get_base <- function(m) {
    if (!identical(base_cache$month, m)) {
      base_cache$base  <<- base_map(m, land)
      base_cache$month <<- m
    }
    base_cache$base
  }
  
  if (PREVIEW_ONLY) {
    t   <- as.POSIXct(PREVIEW_TIME, tz = "UTC")
    dat <- open_data(file_for(t, kind))
    i   <- which.min(abs(as.numeric(dat$times) - as.numeric(t)))
    fields <- dat$get(i)
    dat$close()
    img <- composite(fields, get_base(as.integer(format(t, "%m"))), ocean, ip, luts)
    out <- file.path(LOCAL_DIR, "preview_frame_R.png")
    open_device(out, "png")
    draw_frame(img, stamp_text(dat$times[i]), keys)
    dev.off()
    message("Preview saved: ", out)
    return(invisible())
  }
  
  frame_dir <- file.path(LOCAL_DIR, "frames_R")
  unlink(frame_dir, recursive = TRUE)
  dir.create(frame_dir, recursive = TRUE)
  
  n_written <- 0L
  for (d in seq_along(files)) {
    dat <- open_data(files[d])
    for (i in seq_along(dat$times)) {
      if (!render_this(dat$times[i])) next
      fields <- dat$get(i)
      month  <- as.integer(format(dat$times[i], "%m", tz = "UTC"))
      img    <- composite(fields, get_base(month), ocean, ip, luts)
      n_written <- n_written + 1L
      open_device(file.path(frame_dir, sprintf("frame_%05d.jpg", n_written)), "jpeg")
      draw_frame(img, stamp_text(dat$times[i]), keys)
      dev.off()
    }
    dat$close()
    rm(dat); invisible(gc())
    if (kind == "monthly" || d %% 10 == 0 || d == length(files)) {
      message(sprintf("  file %d/%d  %s  (%d frames)", d, length(files), basename(files[d]), n_written))
    }
  }
  
  frames <- sort(list.files(frame_dir, pattern = "\\.jpg$", full.names = TRUE))
  local_mp4 <- file.path(LOCAL_DIR, OUT_NAME)
  if (file.exists(local_mp4)) file.remove(local_mp4)
  av_encode_video(frames, output = local_mp4, framerate = FPS, codec = "libx264",
                  vfilter = "format=yuv420p", verbose = FALSE)
  
  verify_mp4(local_mp4, n_written)
  unlink(frame_dir, recursive = TRUE)
  
  final_mp4 <- move_to_final(local_mp4)
  message(sprintf("Saved %s (%.0f MB)", final_mp4, file.size(final_mp4) / 1e6))
}

main()