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
# Script: co_surface_platecarree.R
# Description: Generates the surface carbon monoxide animation from NASA MERRA-2
#     hourly CO surface concentration COSC (2006-12-01 to 2006-12-31), rendered
#     in the Plate Carree projection over the December Blue Marble.
# Note: For the data download, see download_merra2_cosc_2006_12.py in the
#     co_surface scripts folder.
#
# Look:
#     Low CO is transparent, so the Blue Marble shows through (navy ocean, white
#     Antarctica). CO fades in as teal/cyan from ~70 ppb and is fully opaque from
#     ~160 ppb, then runs yellow -> red -> black -> grey -> white up to 1600 ppb,
#     matching the original's colorbar (40-1600 ppb, non-linear ticks).
#
# Run:
#     # 1) Set PREVIEW_ONLY <- TRUE, check the PNG
#     # 2) Set PREVIEW_ONLY <- FALSE for the full render
#     caffeinate -i Rscript co_surface_platecarree.R
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
  library(rnaturalearth)
  library(jpeg)
  library(av)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
VARIABLE  <- "CO_surface"                  # ppbv, from download_merra2_cosc_2006_12.py
DRIVE     <- "/Volumes/CLIMATEDATA"
DATA_PATH <- file.path(DRIVE, "co", "merra2_2006_12", "MERRA2_COSC_hourly_2006-12.nc")
BM_DIR    <- file.path(DRIVE, "blue_marble", "2004")
BM_MONTH  <- 12                            # December snow cover, like the original

LOCAL_DIR <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "co_surface", "_render")  # APFS
FINAL_DIR <- file.path(DRIVE, "co", "animations")                                       # exFAT
OUT_NAME  <- "co_surface_dec2006_platecarree_R.mp4"
NE_CACHE  <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "_naturalearth")         # shapefile cache

PREVIEW_ONLY <- FALSE                      # TRUE: save one PNG and stop
PREVIEW_TIME <- "2006-12-07T03:30"         # matches the reference frame (0300)

CENTRAL_LON     <- 165.0                   # Pacific-centered, as in the original
SHOW_COASTLINES <- TRUE
BM_BRIGHTNESS   <- 0.9

# Colorbar ticks from the original; equally spaced on the bar (non-linear scale)
TICKS <- c(40, 100, 160, 200, 400, 600, 800, 1200, 1600)

# (value ppb, color, alpha). Alpha ramps in so clean air shows the Blue Marble.
ANCHORS <- data.frame(
  value = c(40, 70, 100, 140, 160, 175, 190, 200, 300, 400, 600, 800, 1200, 1600),
  color = c("#0b1622", "#12405a", "#1b6f8c", "#3cb4cf", "#7fd6d0", "#f2e03b", "#f5a623",
            "#ec5a1c", "#c0230f", "#1c0705", "#3d3d3d", "#7d7d7d", "#c9c9c9", "#ffffff"),
  alpha = c(0.00, 0.45, 0.75, 0.90, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
  stringsAsFactors = FALSE
)
CMAP_N <- 1024                             # LinearSegmentedColormap N

CB_TITLE    <- "Carbon Monoxide"
CB_SUBTITLE <- "(surface concentration in ppb)"

# Frame layout: 1920 x 1200 px = black top bar, 2:1 map, black bottom bar
W <- 1920
H <- 1200
TOP_PX <- 50
MAP_PX <- 960
DPI <- 100
FPS <- 30                                  # 744 frames -> 24.8 s (original 0:25)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
pt <- function(x) x / 0.75                 # matplotlib points -> R lwd

TICK_POS <- seq(0, 1, length.out = length(TICKS))

# FuncNorm forward: value -> position on the bar (-0.01 / 1.01 outside the range)
fwd <- function(v) {
  p <- approx(TICKS, TICK_POS, xout = v, rule = 2)$y
  p[!is.na(v) & v < TICKS[1]] <- -0.01
  p[!is.na(v) & v > TICKS[length(TICKS)]] <- 1.01
  p
}

# LinearSegmentedColormap equivalent: N RGBA entries sampled at linspace(0, 1, N)
build_cmaps <- function() {
  pos <- fwd(ANCHORS$value)
  rgb_anchor <- t(col2rgb(ANCHORS$color)) / 255
  x <- seq(0, 1, length.out = CMAP_N)
  lut <- cbind(
    r = approx(pos, rgb_anchor[, 1], x)$y,
    g = approx(pos, rgb_anchor[, 2], x)$y,
    b = approx(pos, rgb_anchor[, 3], x)$y,
    a = approx(pos, ANCHORS$alpha, x)$y
  )
  bar <- rgb(lut[, "r"], lut[, "g"], lut[, "b"])  # opaque twin for the colorbar
  list(lut = lut, bar = bar)
}

# Values -> RGBA rows (n x 4). Under (< 40 ppb) and NaN are transparent,
# over (> 1600 ppb) is opaque white, as in the Python cmap.set_under/over/bad.
colorize <- function(vals, lut) {
  p <- fwd(vals) * CMAP_N
  idx <- floor(p)
  idx[!is.na(p) & p == CMAP_N] <- CMAP_N - 1
  out <- matrix(0, length(vals), 4)
  ok <- !is.na(idx) & idx >= 0 & idx < CMAP_N
  out[ok, ] <- lut[idx[ok] + 1, ]
  over <- !is.na(idx) & idx >= CMAP_N
  out[over, ] <- matrix(c(1, 1, 1, 1), sum(over), 4, byrow = TRUE)
  out
}

# Area-average resize of one channel (stand-in for PIL's LANCZOS downsample).
# Weights are banded, so this is a few vector ops per output row/column.
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

load_blue_marble <- function() {
  paths <- sort(list.files(BM_DIR, pattern = sprintf("^world\\.2004%02d.*\\.jpg$", BM_MONTH),
                           full.names = TRUE))
  paths <- paths[!startsWith(basename(paths), "._")]
  if (!length(paths)) stop(sprintf("No Blue Marble image for month %02d in %s", BM_MONTH, BM_DIR))
  message("Loading Blue Marble ", basename(paths[1]), "...")
  img <- readJPEG(paths[1])                                 # [h, w, 3] in 0-1
  if (length(dim(img)) == 2) img <- array(img, c(dim(img), 3))

  # Recenter on CENTRAL_LON, then resize once to the on-screen map size
  w_in <- dim(img)[2]
  shift <- round(((CENTRAL_LON - 180) + 180) / 360 * w_in)
  cols <- ((seq_len(w_in) - 1 + shift) %% w_in) + 1         # np.roll(arr, -shift)
  out <- array(0, c(MAP_PX, W, 3))
  for (k in 1:3) out[, , k] <- resize_channel(img[, cols, k], MAP_PX, W)
  pmin(pmax(out * BM_BRIGHTNESS, 0), 1)
}

# Precompute bilinear lookups from the lat-lon grid onto the map pixels.
# Output pixels run x-fastest, top row first (raster order).
build_resampler <- function(lon, lat) {
  xc <- -180 + (seq_len(W) - 0.5) * 360 / W                 # projected x (deg)
  yc <-   90 - (seq_len(MAP_PX) - 0.5) * 180 / MAP_PX       # latitude, north first
  plon <- rep(xc, times = MAP_PX) + CENTRAL_LON             # geographic lon
  plat <- rep(yc, each = W)

  nlon <- length(lon); nlat <- length(lat)
  dlon <- lon[2] - lon[1]; dlat <- lat[2] - lat[1]

  fi <- ((plon - lon[1]) / dlon) %% nlon                    # wraps across the seam
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

# Alpha-composite the CO layer over the Blue Marble -> one opaque raster
to_raster <- function(vals, bm_rgb, lut) {
  co <- colorize(vals, lut)
  a <- co[, 4]
  rgb_out <- co[, 1:3] * a + bm_rgb * (1 - a)
  m <- matrix(rgb(rgb_out[, 1], rgb_out[, 2], rgb_out[, 3]), nrow = MAP_PX, byrow = TRUE)
  as.raster(m)
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

# sf geometries -> one lon/lat matrix with NA rows separating each line/ring
sf_to_path <- function(geom) {
  xy <- st_coordinates(st_geometry(geom))
  lcols <- grep("^L", colnames(xy))
  grp <- do.call(paste, c(as.data.frame(xy[, lcols, drop = FALSE]), sep = "_"))
  pieces <- split(as.data.frame(xy[, 1:2]), factor(grp, levels = unique(grp)))
  do.call(rbind, lapply(pieces, function(p) rbind(as.matrix(p), c(NA, NA))))
}

# Shift lon/lat paths into the Pacific-centered frame and cut them at the seam,
# adding the crossing point on both map edges (what Cartopy does on its own)
shift_path <- function(lonlat) {
  x <- ((lonlat[, 1] - CENTRAL_LON + 180) %% 360) - 180
  y <- lonlat[, 2]
  dx <- diff(x)
  cut <- which(!is.na(dx) & abs(dx) > 180)
  if (!length(cut)) return(cbind(x, y))

  pieces <- vector("list", 2 * length(cut) + 1)
  prev <- 1
  for (k in seq_along(cut)) {
    j <- cut[k]
    edge <- if (dx[j] < 0) 180 else -180                  # side the line leaves from
    x_next <- x[j + 1] + 2 * edge                          # unwrap the next point
    y_edge <- y[j] + (edge - x[j]) / (x_next - x[j]) * (y[j + 1] - y[j])
    pieces[[2 * k - 1]] <- rbind(cbind(x[prev:j], y[prev:j]), c(edge, y_edge))
    pieces[[2 * k]] <- rbind(c(NA, NA), c(-edge, y_edge))
    prev <- j + 1
  }
  pieces[[2 * length(cut) + 1]] <- cbind(x[prev:length(x)], y[prev:length(y)])
  do.call(rbind, pieces)
}

build_coastlines <- function() {
  if (!SHOW_COASTLINES) return(NULL)
  message("Preparing coastlines...")
  shift_path(sf_to_path(ne_get(50, "coastline", "physical")))
}

stamp_text <- function(t) format(t, "%Y-%m-%d %H%M", tz = "UTC")

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

# Horizontal colorbar block in figure coordinates (cax = [0.71, 0.062, 0.25, 0.026])
draw_colorbar <- function(bar_cols, x0 = 0.71, y0 = 0.062, w = 0.25, h = 0.026) {
  x1 <- x0 + w
  y1 <- y0 + h
  text(0.835, 0.135, CB_TITLE, col = "white", cex = 18 / 12, font = 2, adj = c(0.5, 0.5))
  text(0.835, 0.108, CB_SUBTITLE, col = "#9a9a9a", cex = 1, adj = c(0.5, 0.5))

  rasterImage(as.raster(matrix(bar_cols, nrow = 1)), x0, y0, x1, y1, interpolate = FALSE)

  pad <- 6 / 72 / (H / DPI)                 # 6 pt label pad in figure units
  xt <- x0 + TICK_POS * w
  ends <- c(1, length(TICKS))
  mids <- setdiff(seq_along(TICKS), ends)
  text(xt[mids], y0 - pad, TICKS[mids], col = "#8a8a8a", cex = 1, adj = c(0.5, 1))
  text(xt[ends], y0 - pad, TICKS[ends], col = "white", cex = 15 / 12, adj = c(0.5, 1))
}

draw_frame <- function(img, stamp, coast, bar_cols) {
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), bg = "black", xpd = NA)
  plot.new()

  # Map axes: full width, MAP_PX tall, TOP_PX below the top edge
  map_bottom <- 1 - (TOP_PX + MAP_PX) / H
  par(fig = c(0, 1, map_bottom, map_bottom + MAP_PX / H), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i")
  rasterImage(img, -180, -90, 180, 90, interpolate = FALSE)
  if (!is.null(coast)) {
    par(xpd = FALSE)
    lines(coast, col = "black", lwd = pt(0.45))
    par(xpd = NA)
  }

  # Bottom bar: timestamp (left), colorbar block (right)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(0.115, 0.075, stamp, col = "white", cex = 24 / 12, adj = c(0, 0.5))
  draw_colorbar(bar_cols)
}

render_png <- function(path, img, stamp, coast, bar_cols) {
  png(path, width = W, height = H, res = DPI, bg = "black")
  on.exit(dev.off())
  draw_frame(img, stamp, coast, bar_cols)
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
    f                                      # [lon, lat]
  }

  times <- parse_nc_time(nc)
  n <- length(times)
  message(sprintf("Loaded %s: %d hourly frames, %d x %d", VARIABLE, n, length(lat), length(lon)))

  cm <- build_cmaps()
  bm <- load_blue_marble()
  bm_rgb <- cbind(as.vector(t(bm[, , 1])), as.vector(t(bm[, , 2])),
                  as.vector(t(bm[, , 3])))                  # raster order, x fastest
  rm(bm)
  message("Precomputing map lookup grid...")
  rs <- build_resampler(lon, lat)
  coast <- build_coastlines()

  frame_raster <- function(i) to_raster(resample(rs, read_frame(i)), bm_rgb, cm$lut)

  if (PREVIEW_ONLY) {
    target <- as.POSIXct(PREVIEW_TIME, tz = "UTC", format = "%Y-%m-%dT%H:%M")
    i <- which.min(abs(as.numeric(times) - as.numeric(target)))
    f <- read_frame(i)
    p <- quantile(f, c(1, 50, 90, 99, 99.9) / 100, type = 7, na.rm = TRUE, names = FALSE)
    message(sprintf("%s: p1=%.0f  median=%.0f  p90=%.0f  p99=%.0f  p99.9=%.0f  max=%.0f ppb",
                    stamp_text(times[i]), p[1], p[2], p[3], p[4], p[5], max(f, na.rm = TRUE)))
    out <- file.path(LOCAL_DIR, "preview_frame_R.png")
    render_png(out, frame_raster(i), stamp_text(times[i]), coast, cm$bar)
    message("Preview saved: ", out)
    return(invisible())
  }

  # Render frames to local APFS
  frames_dir <- file.path(LOCAL_DIR, "frames_R")
  unlink(frames_dir, recursive = TRUE)
  dir.create(frames_dir, recursive = TRUE)

  for (i in seq_len(n)) {
    f <- file.path(frames_dir, sprintf("frame_%05d.png", i))
    render_png(f, frame_raster(i), stamp_text(times[i]), coast, cm$bar)
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
