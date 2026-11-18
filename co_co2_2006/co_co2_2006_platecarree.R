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
# Script: co_co2_2006_platecarree.R
# Description: Generates the "CO and CO2, 2006" animation from NASA GMAO MERRA-2
#     column CO (COCL, hourly) and CAMS EGG4 column-mean CO2 (XCO2, 3-hourly)
#     (2006-01-01 to 2006-12-31), rendered in the Plate Carree projection.
#     Column CO2 is drawn in color (377-395 ppmv) and column CO as white plumes
#     (0-6 x 1e18 molecules cm-2) over a base map of blue ocean, grey land and
#     white Antarctica, following NASA GMAO's GEOS-5 Nature Run visualization
#     "Simulated Atmospheric Carbon Concentrations" (SVS 30515).
# Note: R translation of co_co2_2006_platecarree.py; for the data
#     download, see download_co_co2_2006.py in the co_co2_2006 scripts folder.
#
# Frames:
#     One frame per 3-hourly CO2 time step (2,920 frames). The hourly CO fields
#     centered within 45 minutes of that time (MERRA-2 hours are centered at
#     HH:30, so normally the two hours either side) are averaged to match it.
#
# Compositing (back to front):
#     1) Base map: flat blue ocean; land from the January 2004 Blue Marble,
#        desaturated to grey; Antarctica blended to white. Land mask from
#        Natural Earth 50m land (cached as _landmask_R.rds).
#     2) CO2: colormap with opacity ramping up from CO2_ALPHA_LO to
#        CO2_ALPHA_HI, so low-CO2 air lets the base map show.
#     3) CO: white with opacity
#        w = clip((CO - CO_FADE_START) / (6 - CO_FADE_START))^CO_GAMMA.
#
# Run:
#     # 1) Set PREVIEW_ONLY <- TRUE, check the PNGs
#     # 2) Set PREVIEW_ONLY <- FALSE for the full render
#     caffeinate -i Rscript co_co2_2006_platecarree.R
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(jpeg)
  library(av)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
YEAR      <- 2006
DRIVE     <- "/Volumes/CLIMATEDATA"
DATA_DIR  <- file.path(DRIVE, "co_co2", YEAR)                                  # exFAT
WORK_DIR  <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "co_co2_2006")  # APFS
LOCAL_DIR <- file.path(WORK_DIR, "_render")
NE_CACHE  <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "_naturalearth")
MASK_CACHE <- file.path(WORK_DIR, "_landmask_R.rds")

CO_FILE   <- sprintf("MERRA2_COCL_hourly_%d.nc", YEAR)
CO_VAR    <- "CO_column"                   # 1e18 molecules cm-2
CO2_FILE  <- sprintf("CAMS_EGG4_XCO2_3hourly_%d.nc", YEAR)
CO2_VAR   <- "CO2_column"                  # ppmv
BLUE_MARBLE <- file.path(DRIVE, "blue_marble", "2004", "world.200401.3x5400x2700.jpg")

OUT_NAME  <- sprintf("co_co2_%d_platecarree_R.mp4", YEAR)
PREVIEW_ONLY  <- FALSE                     # TRUE: save PNG stills and stop
PREVIEW_TIMES <- c("2006-01-20T00", "2006-04-15T12", "2006-08-20T12", "2006-10-10T12")

# Color scales (the original's colorbars)
CO2_MIN <- 377; CO2_MAX <- 395             # ppmv
CO_MIN  <- 0;   CO_MAX  <- 6               # 1e18 molecules cm-2
CO_FADE_START <- 1.0                       # CO below this (clean air) is invisible
CO_GAMMA      <- 0.8                       # <1 brightens moderate plumes
CO2_ALPHA_LO  <- 379.0                     # ppmv: transparent ...
CO2_ALPHA_HI  <- 382.0                     # ... to fully opaque

# Base map
OCEAN_RGB       <- c(0.11, 0.17, 0.42)     # deep blue, close to the 378-379 color
LAND_BRIGHTNESS <- 0.80                    # scales Blue Marble luminance
LAND_KEEP_COLOR <- 0.15                    # 0 = pure grey, 1 = original colors
ANTARCTIC_WHITE <- 0.85                    # blend of Antarctica toward white

# CO2 colormap (the original's colorbar):
#   black 377 -> navy -> teal -> yellow 383 -> red 385 -> dark red -> purple -> pink 395
CO2_STOPS <- data.frame(
  value = c(377.0, 378.5, 380.0, 381.0, 382.0, 383.0, 384.0,
            385.0, 386.0, 387.0, 388.5, 390.0, 392.0, 395.0),
  color = c("#14142a", "#1c2b6c", "#2a6fb3", "#62c6cf", "#b9e18f", "#f3e15c", "#f4a443",
            "#df4429", "#8f1b1b", "#5b1f5f", "#9a72c2", "#d9b6e9", "#efb5df", "#f8ccee"),
  stringsAsFactors = FALSE
)
CMAP_N <- 1024

# Frame layout: 1920 x 1080 px = full-width 2:1 map over a 120 px black band
W <- 1920
H <- 1080
MAP_PX <- 960
DPI <- 100
FPS <- 30                                  # 2,920 frames -> ~97 s
JPEG_QUALITY <- 95                         # frames as JPEG: 2,920 PNGs would need ~8 GB


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
clip01 <- function(x) pmin(pmax(x, 0), 1)

co2_alpha <- function(v) {
  a <- clip01((v - CO2_ALPHA_LO) / (CO2_ALPHA_HI - CO2_ALPHA_LO))
  a * a * (3 - 2 * a)                      # smoothstep
}

co_weight <- function(v) clip01((v - CO_FADE_START) / (CO_MAX - CO_FADE_START))^CO_GAMMA

# LinearSegmentedColormap equivalent: N RGB rows sampled at linspace(0, 1, N)
build_co2_lut <- function() {
  pos <- (CO2_STOPS$value - CO2_MIN) / (CO2_MAX - CO2_MIN)
  rgb_anchor <- t(col2rgb(CO2_STOPS$color)) / 255
  x <- seq(0, 1, length.out = CMAP_N)
  cbind(approx(pos, rgb_anchor[, 1], x)$y,
        approx(pos, rgb_anchor[, 2], x)$y,
        approx(pos, rgb_anchor[, 3], x)$y)
}

co2_colors <- function(v, lut) {
  idx <- floor(clip01((v - CO2_MIN) / (CO2_MAX - CO2_MIN)) * (CMAP_N - 1)) + 1
  lut[idx, , drop = FALSE]
}

# Area-average resize of one channel (stand-in for PIL's LANCZOS downsample)
area_weights <- function(n_in, n_out) {
  scale <- n_in / n_out
  lapply(seq_len(n_out), function(k) {
    a <- (k - 1) * scale
    b <- k * scale
    i <- seq.int(floor(a) + 1, min(ceiling(b), n_in))
    w <- pmin(i, b) - pmax(i - 1, a)
    list(i = i, w = w / sum(w))
  })
}

resize_channel <- function(ch, out_h, out_w) {
  wy <- area_weights(nrow(ch), out_h)
  tmp <- t(vapply(wy, function(e) colSums(ch[e$i, , drop = FALSE] * e$w),
                  numeric(ncol(ch))))                       # [out_h, in_w]
  wx <- area_weights(ncol(ch), out_w)
  vapply(wx, function(e) as.vector(tmp[, e$i, drop = FALSE] %*% e$w),
         numeric(out_h))                                    # [out_h, out_w]
}

# Land mask at map resolution [MAP_PX, W], TRUE = land, cached after first run
land_mask <- function() {
  if (file.exists(MASK_CACHE)) {
    m <- readRDS(MASK_CACHE)
    if (identical(dim(m), c(MAP_PX, W))) return(m)
  }
  ok <- requireNamespace("sf", quietly = TRUE) &&
        requireNamespace("rnaturalearth", quietly = TRUE)
  land <- if (ok) tryCatch({
    dir.create(NE_CACHE, recursive = TRUE, showWarnings = FALSE)
    tryCatch(                              # reuse the cached copy (.shp or .gpkg)
      rnaturalearth::ne_load(scale = 50, type = "land", category = "physical",
                             destdir = NE_CACHE, returnclass = "sf"),
      error = function(e)
        rnaturalearth::ne_download(scale = 50, type = "land", category = "physical",
                                   destdir = NE_CACHE, load = TRUE, returnclass = "sf"))
  }, error = function(e) { message("  Natural Earth unavailable: ", conditionMessage(e)); NULL })
  if (is.null(land)) return(NULL)

  message("Building land mask from Natural Earth 50m (first run only)...")
  xy <- sf::st_coordinates(sf::st_geometry(land))
  grp <- do.call(paste, c(as.data.frame(xy[, grep("^L", colnames(xy)), drop = FALSE]),
                          sep = "_"))
  rings <- split(as.data.frame(xy[, 1:2]), factor(grp, levels = unique(grp)))
  rings <- rings[vapply(rings, nrow, 1L) >= 3]
  path <- do.call(rbind, lapply(rings, function(p) rbind(as.matrix(p), c(NA, NA))))
  path <- path[-nrow(path), , drop = FALSE]  # polypath rejects a trailing NA row

  # Rasterize by drawing the polygons, white on black, then reading them back
  tmp <- tempfile(fileext = ".jpg")
  jpeg(tmp, width = W, height = MAP_PX, quality = 100, bg = "black", antialias = "none")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i")
  polypath(path[, 1], path[, 2], col = "white", border = NA, rule = "evenodd")
  dev.off()
  m <- readJPEG(tmp)
  if (length(dim(m)) == 3) m <- m[, , 1]
  unlink(tmp)
  m <- m > 0.5
  dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(m, MASK_CACHE)
  m
}

# Base map [MAP_PX * W, 3] in raster order (x fastest, top row first)
build_base <- function() {
  if (!file.exists(BLUE_MARBLE)) stop(BLUE_MARBLE, " not found")
  message("Loading Blue Marble ", basename(BLUE_MARBLE), "...")
  img <- readJPEG(BLUE_MARBLE)
  bm <- array(0, c(MAP_PX, W, 3))
  for (k in 1:3) bm[, , k] <- resize_channel(img[, , k], MAP_PX, W)
  rm(img)
  lum <- 0.299 * bm[, , 1] + 0.587 * bm[, , 2] + 0.114 * bm[, , 3]

  land <- land_mask()
  if (is.null(land)) {                     # Blue Marble ocean is dark blue
    message("  Using a colour-based land mask instead")
    land <- !((bm[, , 3] > bm[, , 1] + 0.03) & (lum < 0.35))
  }
  lat <- 90 - (seq_len(MAP_PX) - 0.5) * 180 / MAP_PX
  ant <- land & matrix(lat < -60, MAP_PX, W)

  out <- matrix(0, MAP_PX * W, 3)
  for (k in 1:3) {
    ch <- lum * LAND_BRIGHTNESS * (1 - LAND_KEEP_COLOR) +
          bm[, , k] * LAND_BRIGHTNESS * LAND_KEEP_COLOR
    ch[ant] <- ch[ant] * (1 - ANTARCTIC_WHITE) + ANTARCTIC_WHITE
    ch[!land] <- OCEAN_RGB[k]
    out[, k] <- as.vector(t(clip01(ch)))
  }
  out
}

# Precompute bilinear lookups from a lat-lon grid onto the map pixels.
# Output pixels run x-fastest, top row first (raster order).
build_resampler <- function(lon, lat) {
  plon <- rep(-180 + (seq_len(W) - 0.5) * 360 / W, times = MAP_PX)
  plat <- rep(90 - (seq_len(MAP_PX) - 0.5) * 180 / MAP_PX, each = W)

  nlon <- length(lon); nlat <- length(lat)
  dlon <- lon[2] - lon[1]; dlat <- lat[2] - lat[1]

  fi <- ((plon - lon[1]) / dlon) %% nlon                    # wraps across the dateline
  i0 <- floor(fi)
  wx <- fi - i0
  i0 <- i0 %% nlon
  i1 <- (i0 + 1) %% nlon

  fj <- pmin(pmax((plat - lat[1]) / dlat, 0), nlat - 1)     # edge rows clamp
  j0 <- pmin(floor(fj), nlat - 2)
  wy <- fj - j0
  j1 <- j0 + 1

  list(i00 = j0 * nlon + i0 + 1, i01 = j0 * nlon + i1 + 1,
       i10 = j1 * nlon + i0 + 1, i11 = j1 * nlon + i1 + 1,
       wx = wx, wy = wy)
}

resample <- function(rs, field) {
  f <- as.vector(field)                    # [lon, lat] layout, lon fastest
  top <- (1 - rs$wx) * f[rs$i00] + rs$wx * f[rs$i01]
  bot <- (1 - rs$wx) * f[rs$i10] + rs$wx * f[rs$i11]
  (1 - rs$wy) * top + rs$wy * bot
}

# Base <- CO2 (alpha) <- CO (white, alpha) -> one opaque raster
to_raster <- function(co2, co, base, lut) {
  a2 <- co2_alpha(co2)
  w <- co_weight(co)
  rgb2 <- co2_colors(co2, lut)
  mix <- rgb2 * a2 + base * (1 - a2)
  out <- clip01(w + mix * (1 - w))
  as.raster(matrix(rgb(out[, 1], out[, 2], out[, 3]), nrow = MAP_PX, byrow = TRUE))
}

parse_nc_time <- function(nc) {
  units <- ncatt_get(nc, "time", "units")$value
  m <- regmatches(units, regexec("^\\s*(\\w+) since (.+?)\\s*$", units))[[1]]
  if (length(m) != 3) stop("Cannot parse time units: ", units)
  scale <- switch(tolower(m[2]), seconds = 1, minutes = 60, hours = 3600,
                  days = 86400, stop("Unknown time unit: ", m[2]))
  origin <- as.POSIXct(m[3], tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS",
                                      "%Y-%m-%d %H:%M", "%Y-%m-%d"))
  origin + as.numeric(ncvar_get(nc, "time")) * scale
}

# Opens one variable and serves it one day-sized block at a time
open_series <- function(path, var) {
  nc <- nc_open(path)
  dims <- vapply(nc$var[[var]]$dim, function(d) d$name, "")
  if (!identical(dims, c("lon", "lat", "time"))) {
    stop(basename(path), ": unexpected dimension order ", paste(dims, collapse = ", "))
  }
  lon <- as.numeric(ncvar_get(nc, "lon"))
  lat <- as.numeric(ncvar_get(nc, "lat"))
  if (any(diff(lon) <= 0) || any(diff(lat) <= 0)) {
    stop(basename(path), ": lon/lat must be ascending (as written by the download script)")
  }
  times <- parse_nc_time(nc)
  per_day <- round(86400 / as.numeric(difftime(times[2], times[1], units = "secs")))
  cache <- new.env()
  get <- function(i) {                     # i is 1-based
    day <- (i - 1) %/% per_day
    key <- as.character(day)
    if (is.null(cache[[key]])) {
      keys <- ls(cache)
      if (length(keys) >= 2) rm(list = keys[which.min(as.numeric(keys))], envir = cache)
      s <- day * per_day + 1
      cnt <- min(per_day, length(times) - s + 1)
      cache[[key]] <- ncvar_get(nc, var, start = c(1, 1, s), count = c(-1, -1, cnt),
                                collapse_degen = FALSE)
    }
    cache[[key]][, , i - day * per_day]
  }
  list(nc = nc, lon = lon, lat = lat, times = times, per_day = per_day, get = get)
}

locate <- function(name) {
  for (d in c(DATA_DIR, WORK_DIR)) {
    p <- file.path(d, name)
    if (file.exists(p)) return(p)
  }
  stop("Cannot find ", name, " in ", DATA_DIR, " or ", WORK_DIR)
}

# Horizontal colorbar in figure coordinates, label above, ticks below
draw_colorbar <- function(cols, x0, w, ticks, labels, title) {
  band <- (H - MAP_PX) / H
  y0 <- band * 0.38
  y1 <- y0 + band * 0.20
  rasterImage(as.raster(matrix(cols, nrow = 1)), x0, y0, x0 + w, y1, interpolate = TRUE)
  text(x0, band * 0.80, title, col = "white", cex = 11 / 12, adj = c(0, 0.5))
  xt <- x0 + (ticks - ticks[1]) / (ticks[length(ticks)] - ticks[1]) * w
  text(xt, y0 - band * 0.06, labels, col = "white", cex = 9 / 12, adj = c(0.5, 1))
}

draw_frame <- function(img, t, bars) {
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), bg = "black", xpd = NA)
  plot.new()

  # Map: full width, top MAP_PX rows
  band <- (H - MAP_PX) / H
  par(fig = c(0, 1, band, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i")
  rasterImage(img, -180, -90, 180, 90, interpolate = FALSE)

  # Bottom band: date block (left), colorbars (right)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  date_txt <- format(t, "%Y / %m / %d", tz = "UTC")
  text(0.125, band * 0.60, date_txt, col = "white", cex = 40 / 12, font = 2, adj = c(0, 0.5))
  x_hour <- 0.125 + bars$date_w + 0.008
  text(x_hour, band * 0.60, format(t, "%Hz", tz = "UTC"), col = "#bbbbbb",
       cex = 20 / 12, adj = c(0, 0.5))
  text(0.125, band * 0.15, "SDSU Climate Informatics Lab", col = "white",
       cex = 13 / 12, font = 2, adj = c(0, 0.5))

  draw_colorbar(bars$co, 0.385, 0.29, seq(CO_MIN, CO_MAX, by = 0.6),
                sprintf("%.1f", seq(CO_MIN, CO_MAX, by = 0.6)),
                expression("Carbon Monoxide Column Abundance [10"^18 * " molec cm"^-2 * "]"))
  draw_colorbar(bars$co2, 0.695, 0.29, seq(CO2_MIN, CO2_MAX, by = 2),
                sprintf("%.0f", seq(CO2_MIN, CO2_MAX, by = 2)),
                "Carbon Dioxide Column Concentration [ppmv]")
}

# Widest date of the year, so the hour label never overlaps it
date_width <- function(times) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = W, height = H, res = DPI)
  on.exit({ dev.off(); unlink(tmp) })
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  d <- unique(format(times, "%Y / %m / %d", tz = "UTC"))
  max(strwidth(d, cex = 40 / 12, font = 2))
}

build_bars <- function(lut, times) {
  v <- seq(CO2_MIN, CO2_MAX, length.out = 256)
  co2_rgb <- co2_colors(v, lut) * co2_alpha(v)             # as seen over black
  g <- co_weight(seq(CO_MIN, CO_MAX, length.out = 256))
  list(co2 = rgb(co2_rgb[, 1], co2_rgb[, 2], co2_rgb[, 3]),
       co = rgb(g, g, g),
       date_w = date_width(times))
}

render_frame <- function(path, img, t, bars) {
  if (grepl("\\.png$", path)) {
    png(path, width = W, height = H, res = DPI, bg = "black")
  } else {
    jpeg(path, width = W, height = H, res = DPI, bg = "black", quality = JPEG_QUALITY)
  }
  on.exit(dev.off())
  draw_frame(img, t, bars)
}

verify_mp4 <- function(path, expected_frames) {
  if (!file.exists(path) || file.size(path) == 0) stop(basename(path), " is missing or empty")
  info <- av_media_info(path)
  n_frames <- info$video$frames
  # av_media_info can over-report by one frame (container metadata), so allow +/-1
  if (is.null(n_frames) || is.na(n_frames) || abs(n_frames - expected_frames) > 1) {
    stop(sprintf("%s: %s frames, expected %d", basename(path), n_frames, expected_frames))
  }
  message(sprintf("Verified %s: ~%d frames, %.1f s (expected %.1f s)", basename(path),
                  n_frames, info$duration, expected_frames / FPS))
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  dir.create(LOCAL_DIR, recursive = TRUE, showWarnings = FALSE)
  co <- open_series(locate(CO_FILE), CO_VAR)
  co2 <- open_series(locate(CO2_FILE), CO2_VAR)
  on.exit({ nc_close(co$nc); nc_close(co2$nc) }, add = TRUE)

  fmt <- function(t) format(t, "%Y-%m-%d %H:%M", tz = "UTC")
  message(sprintf("CO : %d hourly steps, %s to %s", length(co$times),
                  fmt(co$times[1]), fmt(co$times[length(co$times)])))
  message(sprintf("CO2: %d 3-hourly steps, %s to %s", length(co2$times),
                  fmt(co2$times[1]), fmt(co2$times[length(co2$times)])))

  # Match each CO2 time to the CO hours centered within 45 minutes of it
  co_s <- as.numeric(co$times)
  matches <- lapply(as.numeric(co2$times), function(t) which(abs(co_s - t) <= 45 * 60))
  if (any(lengths(matches) == 0)) {
    stop("No CO hour within 45 min of ", fmt(co2$times[which(lengths(matches) == 0)[1]]))
  }
  n <- length(matches)
  message(sprintf("Matched %d frames (%d averaged from 2 CO hours)", n, sum(lengths(matches) == 2)))

  # Scale check against the original's colorbars (one field per week)
  wk <- seq(1, n, by = co2$per_day * 7)
  s <- unlist(lapply(wk, co2$get))
  p <- quantile(s, c(0.01, 0.5, 0.99), names = FALSE)
  message(sprintf("CO2 weekly sample: p1=%.1f median=%.1f p99=%.1f ppmv (scale %.0f-%.0f)",
                  p[1], p[2], p[3], CO2_MIN, CO2_MAX))

  lut <- build_co2_lut()
  base <- build_base()
  message("Precomputing map lookup grids...")
  rs_co <- build_resampler(co$lon, co$lat)
  rs_co2 <- build_resampler(co2$lon, co2$lat)
  bars <- build_bars(lut, co2$times)

  frame_raster <- function(k) {
    co_field <- Reduce(`+`, lapply(matches[[k]], co$get)) / length(matches[[k]])
    to_raster(resample(rs_co2, co2$get(k)), resample(rs_co, co_field), base, lut)
  }

  if (PREVIEW_ONLY) {
    for (pt in PREVIEW_TIMES) {
      target <- as.POSIXct(pt, tz = "UTC", format = "%Y-%m-%dT%H")
      k <- which.min(abs(as.numeric(co2$times) - as.numeric(target)))
      out <- file.path(LOCAL_DIR, format(co2$times[k], "preview_R_%Y%m%d_%H.png", tz = "UTC"))
      render_frame(out, frame_raster(k), co2$times[k], bars)
      message("  wrote ", out)
    }
    return(invisible())
  }

  # Render frames to local APFS
  frames_dir <- file.path(LOCAL_DIR, "frames_R")
  unlink(frames_dir, recursive = TRUE)
  dir.create(frames_dir, recursive = TRUE)
  t0 <- Sys.time()
  for (k in seq_len(n)) {
    render_frame(file.path(frames_dir, sprintf("frame_%05d.jpg", k)),
                 frame_raster(k), co2$times[k], bars)
    if ((k - 1) %% 40 == 0 || k == n) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf("  frame %d/%d  %s  (%.1f min elapsed, ~%.0f min left)",
                      k, n, fmt(co2$times[k]), el, el / k * (n - k)))
    }
  }

  # Encode
  local_mp4 <- file.path(LOCAL_DIR, OUT_NAME)
  if (file.exists(local_mp4)) file.remove(local_mp4)
  frames <- list.files(frames_dir, pattern = "^frame_\\d+\\.jpg$", full.names = TRUE)
  av_encode_video(frames, output = local_mp4, framerate = FPS,
                  codec = "libx264", vfilter = "format=yuv420p", verbose = FALSE)
  verify_mp4(local_mp4, n)
  unlink(frames_dir, recursive = TRUE)

  # Move to exFAT (copy + remove: file.rename fails across filesystems)
  if (!dir.exists(DRIVE)) {
    message(DRIVE, " not mounted - video left at ", local_mp4)
    return(invisible())
  }
  dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
  final_mp4 <- file.path(DATA_DIR, OUT_NAME)
  if (!file.copy(local_mp4, final_mp4, overwrite = TRUE) ||
      file.size(final_mp4) != file.size(local_mp4)) {
    stop("Copy to ", DATA_DIR, " failed; local file kept at ", local_mp4)
  }
  file.remove(local_mp4)
  message("Saved ", final_mp4)
}

if (sys.nframe() == 0L) main()
