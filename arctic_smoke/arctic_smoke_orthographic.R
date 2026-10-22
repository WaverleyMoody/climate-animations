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
# Script: arctic_smoke_orthographic.R
# Description: Generates the Arctic boreal forest wildfire smoke animation from
#     NASA MERRA-2 hourly total aerosol optical depth TOTEXTTAU (2019-06-01 to
#     2019-08-31), rendered in the Orthographic projection.
# Note: For the data download, see download_merra2_aod_2019_jja.py in the
#     arctic_smoke scripts folder.
#
# Run:
#     # 1) Set PREVIEW_ONLY <- TRUE, check the PNG
#     # 2) Set PREVIEW_ONLY <- FALSE for the full render
#     caffeinate -i Rscript arctic_smoke_orthographic.R
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
  library(rnaturalearth)
  library(av)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
VARIABLE  <- "TOTEXTTAU"
DATA_PATH <- "/Volumes/CLIMATEDATA/smoke/merra2_2019_jja/MERRA2_TOTEXTTAU_hourly_2019-06-08.nc"

LOCAL_DIR <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "arctic_smoke", "_render")  # local APFS
FINAL_DIR <- "/Volumes/CLIMATEDATA/smoke/animations"                                       # exFAT
OUT_NAME  <- "arctic_smoke_jja2019_orthographic_R.mp4"
NE_CACHE  <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "_naturalearth")           # shapefile cache

PREVIEW_ONLY <- FALSE                      # TRUE: save one PNG and stop
PREVIEW_TIME <- "2019-06-13T03:30"         # matches the reference frame (04z)

# Globe view, measured from the reference frame
CENTRAL_LON <- -115.0
CENTRAL_LAT <- 70.0
LAT_CIRCLES <- c(63, 80)                   # dashed parallels in the reference
EARTH_R     <- 6378137                     # sphere radius (m) for the projection

# Resolution of the precomputed projection grid (pixels across the globe)
GRID_N <- 1400

# Colorbar levels from the original (log-spaced, extended both ends)
LEVELS <- c(0.10, 0.15, 0.21, 0.31, 0.45, 0.66, 0.97, 1.40, 2.10, 3.00)

BG_COLOR     <- "#fdf5e8"                  # cream globe background
CMAP_ANCHORS <- c("#fbe9cf", "#f6cf9c", "#eeac72", "#e2864f",
                  "#cd5f34", "#ad3c21", "#842312", "#5a0e06")
OVER_COLOR   <- "#3d0703"

# Line styles: light coastlines worldwide, dark outlines for US/Canada + states
LIGHT_LINE     <- "#9a9a9a"
DARK_LINE      <- "#3a3a3a"
DARK_COUNTRIES <- c("USA", "CAN")

SHOW_US_CITIES <- TRUE
CITY_MIN_POP   <- 50000

TITLE      <- "Smoke from Arctic Boreal Forest Wildfires, June\u2013August 2019"
SUBTITLE   <- "NASA MERRA-2 hourly aerosol optical depth (550 nm)"
CBAR_LABEL <- "Dust optical depth (550nm)"

FPS     <- 35                              # 2208 frames -> 63 s, same as original
FIGSIZE <- c(12, 10.4)                     # 1800 x 1560 px at DPI 150
DPI     <- 150

ORTHO_CRS <- sprintf("+proj=ortho +lat_0=%g +lon_0=%g +R=%g +units=m +no_defs",
                     CENTRAL_LAT, CENTRAL_LON, EARTH_R)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
pt <- function(x) x / 0.75                 # matplotlib points -> R lwd

build_palette <- function() {
  # Same as sampling the smooth colormap at linspace(0, 1, n_bins)
  colorRampPalette(CMAP_ANCHORS)(length(LEVELS) - 1)
}

ne_get <- function(scale, type, category) {
  dir.create(NE_CACHE, recursive = TRUE, showWarnings = FALSE)
  shp <- file.path(NE_CACHE, sprintf("ne_%dm_%s.shp", scale, type))
  if (file.exists(shp)) {
    ne_load(scale = scale, type = type, category = category,
            destdir = NE_CACHE, returnclass = "sf")
  } else {
    ne_download(scale = scale, type = type, category = category,
                destdir = NE_CACHE, load = TRUE, returnclass = "sf")
  }
}

attr_col <- function(x, name) {
  nm <- names(x)[toupper(names(x)) == name]
  if (!length(nm)) stop("Column ", name, " not found")
  x[[nm[1]]]
}

# sf geometries -> one lon/lat matrix with NA rows separating each line/ring
sf_to_path <- function(geom) {
  xy <- st_coordinates(st_geometry(geom))
  lcols <- grep("^L", colnames(xy))
  grp <- do.call(paste, c(as.data.frame(xy[, lcols, drop = FALSE]), sep = "_"))
  pieces <- split(as.data.frame(xy[, 1:2]), factor(grp, levels = unique(grp)))
  do.call(rbind, lapply(pieces, function(p) rbind(as.matrix(p), c(NA, NA))))
}

# Forward-project lon/lat; points on the far side of the globe become NA,
# which breaks the line there (same effect as Cartopy's clipping)
project_path <- function(lonlat, crs = ORTHO_CRS) {
  out <- matrix(NA_real_, nrow(lonlat), 2)
  ok <- !is.na(lonlat[, 1]) & !is.na(lonlat[, 2])
  p <- sf_project("EPSG:4326", crs, lonlat[ok, , drop = FALSE],
                  keep = TRUE, warn = FALSE)
  p[!is.finite(p)] <- NA
  out[ok, ] <- p
  out
}

us_city_points <- function() {
  tryCatch({
    pp <- ne_get(10, "populated_places", "cultural")
    pop <- attr_col(pp, "POP_MAX")
    keep <- attr_col(pp, "ADM0_A3") == "USA" & !is.na(pop) & pop >= CITY_MIN_POP
    xy <- project_path(st_coordinates(st_geometry(pp)[keep])[, 1:2, drop = FALSE])
    xy[complete.cases(xy), , drop = FALSE]
  }, error = function(e) {
    message("City layer skipped (", conditionMessage(e), ")")
    matrix(numeric(0), ncol = 2)
  })
}

build_layers <- function() {
  message("Preparing map layers...")
  coast     <- ne_get(50, "coastline", "physical")
  lakes     <- ne_get(50, "lakes", "physical")
  countries <- ne_get(50, "admin_0_countries", "cultural")
  countries <- countries[attr_col(countries, "ADM0_A3") %in% DARK_COUNTRIES, ]
  states    <- ne_get(50, "admin_1_states_provinces_lakes", "cultural")
  states    <- states[attr_col(states, "ADM0_A3") == "USA", ]
  
  circle_lons <- seq(-180, 180, length.out = 721)
  circles <- do.call(rbind, lapply(LAT_CIRCLES, function(lt)
    rbind(cbind(circle_lons, lt), c(NA, NA))))
  
  list(
    coast     = project_path(sf_to_path(coast)),
    lakes     = project_path(sf_to_path(lakes)),
    countries = project_path(sf_to_path(countries)),
    states    = project_path(sf_to_path(states)),
    circles   = project_path(circles),
    cities    = if (SHOW_US_CITIES) us_city_points() else matrix(numeric(0), ncol = 2)
  )
}

stamp_text <- function(t) format(t, "%Y-%m-%d %H:%Mz", tz = "UTC")

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

# Precompute bilinear lookups from the lat-lon grid onto the globe image.
# The projection math runs once, so each frame is just an index gather.
# Grid vectors run x-fastest, y from the bottom (like origin="lower").
build_resampler <- function(lon, lat, n, r = EARTH_R) {
  centers <- -r + (seq_len(n) - 0.5) * (2 * r / n)
  xx <- rep(centers, times = n)
  yy <- rep(centers, each = n)
  inside <- (xx^2 + yy^2) < r^2
  
  ll <- sf_project(ORTHO_CRS, "EPSG:4326", cbind(xx[inside], yy[inside]),
                   keep = TRUE, warn = FALSE)
  plon <- ll[, 1]
  plat <- ll[, 2]
  
  nlon <- length(lon); nlat <- length(lat)
  dlon <- lon[2] - lon[1]; dlat <- lat[2] - lat[1]
  
  fi <- ((plon - lon[1]) / dlon) %% nlon
  i0 <- floor(fi)
  wx <- fi - i0
  i0 <- i0 %% nlon
  i1 <- (i0 + 1) %% nlon
  
  fj <- (plat - lat[1]) / dlat
  valid <- is.finite(fi) & is.finite(fj) & fj >= 0 & fj <= nlat - 1
  j0 <- pmin(pmax(floor(fj), 0), nlat - 2)
  wy <- pmin(pmax(fj - j0, 0), 1)
  j1 <- j0 + 1
  
  idx <- function(j, i) (j * nlon + i + 1)[valid]   # 1-based, lon varies fastest
  
  list(n = n, inside = inside, target = which(inside)[valid],
       i00 = idx(j0, i0), i01 = idx(j0, i1),
       i10 = idx(j1, i0), i11 = idx(j1, i1),
       wx = wx[valid], wy = wy[valid])
}

resample <- function(rs, field) {
  f <- as.vector(field)                    # [lon, lat] matrix, lon fastest
  top <- (1 - rs$wx) * f[rs$i00] + rs$wx * f[rs$i01]
  bot <- (1 - rs$wx) * f[rs$i10] + rs$wx * f[rs$i11]
  out <- rep(NA_real_, rs$n^2)
  out[rs$target] <- (1 - rs$wy) * top + rs$wy * bot
  out
}

# Values -> color raster (BoundaryNorm equivalent: under = BG, >= 3.00 = over)
to_raster <- function(rs, vals, pal) {
  lut <- c(BG_COLOR, pal, OVER_COLOR)
  cols <- lut[findInterval(vals, LEVELS) + 1]
  cols[is.na(vals) & rs$inside] <- BG_COLOR
  cols[!rs$inside] <- "transparent"
  m <- t(matrix(cols, nrow = rs$n))        # rows = y (bottom -> top)
  as.raster(m[rs$n:1, ])                   # raster rows run top -> bottom
}

draw_colorbar <- function(pal, x0 = 0.895, y0 = 0.27, w = 0.018, h = 0.46) {
  x1 <- x0 + w
  y1 <- y0 + h
  nb <- length(pal)
  edges <- seq(y0, y1, length.out = nb + 1)
  ext <- 0.05 * h                          # matplotlib default extendfrac
  
  rect(x0, edges[-(nb + 1)], x1, edges[-1], col = pal, border = NA)
  polygon(c(x0, x1, (x0 + x1) / 2), c(y0, y0, y0 - ext), col = BG_COLOR, border = NA)
  polygon(c(x0, x1, (x0 + x1) / 2), c(y1, y1, y1 + ext), col = OVER_COLOR, border = NA)
  polygon(c(x0, (x0 + x1) / 2, x1, x1, (x0 + x1) / 2, x0),
          c(y0, y0 - ext, y0, y1, y1 + ext, y1), col = NA, border = "black", lwd = pt(0.8))
  
  tick <- 3.5 / 72 / FIGSIZE[1]            # 3.5 pt tick length in figure units
  segments(x1, edges, x1 + tick, edges, lwd = pt(0.8))
  text(x1 + 2 * tick, edges, sprintf("%.2f", LEVELS), adj = c(0, 0.5), cex = 13 / 12)
  
  pad <- 8 / 72 / FIGSIZE[1]               # labelpad = 8 pt
  text(x0 - pad, (y0 + y1) / 2, CBAR_LABEL, srt = 90, adj = c(0.5, 0), cex = 14 / 12)
}

draw_frame <- function(img, stamp, layers, pal, r = EARTH_R) {
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), bg = "white", xpd = NA)
  plot.new()
  
  # Globe axes: [0.01, 0.03, 0.85, 0.86] in figure coordinates
  par(fig = c(0.01, 0.86, 0.03, 0.89), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-r, r), ylim = c(-r, r), asp = 1, xaxs = "i", yaxs = "i")
  
  th <- seq(0, 2 * pi, length.out = 721)
  polygon(r * cos(th), r * sin(th), col = BG_COLOR, border = NA)
  rasterImage(img, -r, -r, r, r, interpolate = FALSE)
  
  # Light layers worldwide
  lines(layers$coast, col = LIGHT_LINE, lwd = pt(0.45))
  lines(layers$lakes, col = LIGHT_LINE, lwd = pt(0.35))
  
  # Dark outlines for the US and Canada, plus US states
  lines(layers$countries, col = DARK_LINE, lwd = pt(0.6))
  lines(layers$states, col = DARK_LINE, lwd = pt(0.45))
  
  # Dashed latitude circles
  lines(layers$circles, col = "#777777", lwd = pt(0.7), lty = "44")
  
  if (nrow(layers$cities)) {
    points(layers$cities, pch = 16, cex = 0.15, col = "black")
  }
  
  # Figure-level text and colorbar
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(0.435, 0.955, TITLE, adj = c(0.5, 0), cex = 17 / 12, font = 2)
  text(0.435, 0.925, SUBTITLE, adj = c(0.5, 0), cex = 1, col = "#444444")
  draw_colorbar(pal)
  text(0.93, 0.03, stamp, adj = c(1, 0), cex = 13 / 12)
}

render_png <- function(path, img, stamp, layers, pal) {
  png(path, width = FIGSIZE[1] * DPI, height = FIGSIZE[2] * DPI, res = DPI)
  on.exit(dev.off())
  draw_frame(img, stamp, layers, pal)
}

verify_mp4 <- function(path, expected_frames) {
  if (!file.exists(path) || file.size(path) == 0) stop(basename(path), " is missing or empty")
  info <- av_media_info(path)
  n_frames <- info$video$frames
  if (is.null(n_frames) || is.na(n_frames) || n_frames != expected_frames) {
    stop(sprintf("%s: %s frames, expected %d", basename(path), n_frames, expected_frames))
  }
  message(sprintf("Verified %s: %d frames, %.1f s", basename(path), n_frames, info$duration))
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  if (!file.exists(DATA_PATH)) stop(DATA_PATH, " not found. Is /Volumes/CLIMATEDATA mounted?")
  dir.create(LOCAL_DIR, recursive = TRUE, showWarnings = FALSE)
  
  nc <- nc_open(DATA_PATH)
  on.exit(nc_close(nc), add = TRUE)
  
  dims <- vapply(nc$var[[VARIABLE]]$dim, function(d) d$name, "")
  if (!identical(dims[1:2], c("lon", "lat"))) {
    stop("Unexpected dimension order: ", paste(dims, collapse = ", "))
  }
  
  lon <- as.numeric(ncvar_get(nc, "lon"))
  lat <- as.numeric(ncvar_get(nc, "lat"))
  flip_lat <- lat[1] > lat[length(lat)]
  if (flip_lat) lat <- rev(lat)
  if (any(diff(lon) <= 0)) stop("Longitude axis is not ascending; check the input file")
  
  read_frame <- function(i) {
    f <- ncvar_get(nc, VARIABLE, start = c(1, 1, i), count = c(-1, -1, 1))
    if (flip_lat) f <- f[, ncol(f):1]
    f
  }
  
  times <- parse_nc_time(nc)
  n <- length(times)
  message(sprintf("Loaded %s: %d hourly frames, %d x %d", VARIABLE, n, length(lat), length(lon)))
  
  message("Precomputing orthographic lookup grid...")
  rs <- build_resampler(lon, lat, GRID_N)
  layers <- build_layers()
  pal <- build_palette()
  
  if (PREVIEW_ONLY) {
    target <- as.POSIXct(PREVIEW_TIME, tz = "UTC", format = "%Y-%m-%dT%H:%M")
    i <- which.min(abs(as.numeric(times) - as.numeric(target)))
    out <- file.path(LOCAL_DIR, "preview_frame_R.png")
    render_png(out, to_raster(rs, resample(rs, read_frame(i)), pal),
               stamp_text(times[i]), layers, pal)
    message("Preview saved: ", out)
    return(invisible())
  }
  
  # Render frames to local APFS
  frames_dir <- file.path(LOCAL_DIR, "frames_R")
  unlink(frames_dir, recursive = TRUE)
  dir.create(frames_dir, recursive = TRUE)
  
  for (i in seq_len(n)) {
    f <- file.path(frames_dir, sprintf("frame_%05d.png", i))
    render_png(f, to_raster(rs, resample(rs, read_frame(i)), pal),
               stamp_text(times[i]), layers, pal)
    if ((i - 1) %% 24 == 0) message(sprintf("  frame %d/%d  %s", i, n, stamp_text(times[i])))
  }
  
  # Encode
  local_mp4 <- file.path(LOCAL_DIR, OUT_NAME)
  if (file.exists(local_mp4)) file.remove(local_mp4)
  frames <- list.files(frames_dir, pattern = "^frame_\\d+\\.png$", full.names = TRUE)
  av_encode_video(frames, output = local_mp4, framerate = FPS,
                  codec = "libx264", vfilter = "format=yuv420p", verbose = FALSE)
  
  verify_mp4(local_mp4, n)
  unlink(frames_dir, recursive = TRUE)
  
  # Move to exFAT (copy + remove: file.rename fails across filesystems)
  dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
  final_mp4 <- file.path(FINAL_DIR, OUT_NAME)
  if (!file.copy(local_mp4, final_mp4, overwrite = TRUE) ||
      file.size(final_mp4) != file.size(local_mp4)) {
    stop("Copy to ", FINAL_DIR, " failed; local file kept at ", local_mp4)
  }
  file.remove(local_mp4)
  message("Saved ", final_mp4)
}

if (sys.nframe() == 0L) main()