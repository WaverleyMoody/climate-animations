# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: animate_snow_sea_ice.R
#
# Description: Generates the snow cover and sea ice animation from IMS 4 km NH snow/ice
#              (NSIDC G02156) and NOAA/NSIDC Sea Ice Concentration CDR V6 south (G02202_V6)
#              (2019-2021), rendered as side-by-side north and south polar orthographic globes
#              over Blue Marble Next Generation (2004) imagery. R translation of
#              animate_snow_sea_ice.py.
#
# Note: Adapted from NASA SVS 4995 ("Global Snow Cover and Sea Ice Cycle at Both Poles").
#
# Approach: every output pixel of each globe is mapped once (at startup) to its source pixel in
# the Blue Marble image, the IMS grid, and the G02202 grid. Each frame is then pure vector
# indexing, so memory stays low and rendering is fast.
#
# Packages: ncdf4, sf, jpeg, av
#           install.packages(c("ncdf4", "sf", "jpeg", "av"))
#
# Note: Frames are written as PNGs to local disk and encoded with av at the end (~600 MB for
#       the default 5-day step; deleted after the video is verified and moved).
#       Run with: caffeinate -i Rscript animate_snow_sea_ice.R

suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
  library(jpeg)
  library(av)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
START <- as.Date("2019-01-01")
END   <- as.Date("2021-12-31")
FRAME_STEP_DAYS <- 5           # 5 -> ~220 frames; 1 -> daily (1,096 frames)
FPS <- 6                       # 6 fps at 5-day steps ~ 37 s; for daily frames use ~18
PREVIEW_ONLY <- FALSE          # TRUE = write one PNG and stop
PREVIEW_DATE <- as.Date("2019-08-05")

DRIVE   <- "/Volumes/CLIMATEDATA"
IMS_DIR <- file.path(DRIVE, "ims_4km")
SIC_DIR <- file.path(DRIVE, "sea_ice_concentration_south")
BM_DIR  <- file.path(DRIVE, "blue_marble", "2004")

SCRIPT_DIR <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "snow_and_sea_ice")
OUT_NAME   <- "snow_sea_ice_2019_2021_R.mp4"
LOCAL_OUT  <- file.path(SCRIPT_DIR, OUT_NAME)                       # APFS
FINAL_OUT  <- file.path(DRIVE, "snow_and_sea_ice", OUT_NAME)        # exFAT
FRAME_DIR  <- file.path(SCRIPT_DIR, "_frames_tmp_R")                # APFS scratch

# Blue Marble month per globe: each hemisphere's snow-minimum month, so all snow
# and ice on screen comes from the 2019-2021 data, not from 2004 imagery.
NORTH_BM_MONTH <- 8
SOUTH_BM_MONTH <- 2

# Frame layout (1920x1080)
W <- 1920; H <- 1080
R <- 420                                                   # globe radius, px
CENTERS <- list(north = c(510, 555), south = c(1410, 555))
DATE_XY <- c(1880, 1045)                                   # lower-right date label (right-aligned)
FONT_FAMILY <- "Helvetica Neue"   # cairo falls back to the default sans font if missing

# Colours / blending
SNOW_RGB <- c(248, 248, 252)
ICE_RGB  <- c(232, 240, 250)
SNOW_ALPHA <- 0.90
ICE_ALPHA  <- 0.95
SIC_FULL   <- 60.0          # concentration (%) at which southern ice is fully opaque
SIC_GAMMA  <- 0.7
GRID_RGB   <- c(45, 45, 50)
GRID_ALPHA <- 0.35
HAZE_RGB   <- c(170, 195, 230)
HAZE_POWER <- 8; HAZE_STRENGTH <- 0.55

IMS_VAR <- "IMS_Surface_Values"     # 1 water, 2 snow-free land, 3 sea/lake ice, 4 snow on land
SIC_VAR <- "cdr_seaice_conc"        # uint8 0-100 %, 255 = land / fill


# ---------------------------------------------------------------------------
# File indexing
# ---------------------------------------------------------------------------
# Named character vector of paths, keyed by integer day (as.integer(Date)).
index_files <- function(root, regex, fmt) {
  paths <- list.files(root, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)
  paths <- paths[!startsWith(basename(paths), "._")]
  m <- regmatches(basename(paths), regexec(regex, basename(paths)))
  hit <- lengths(m) == 2
  if (!any(hit)) stop("No files found under ", root)
  days <- as.integer(as.Date(vapply(m[hit], `[`, "", 2), format = fmt))
  setNames(paths[hit], days)
}

nearest_file <- function(files, d, max_days = 3) {
  d <- as.integer(d)
  for (k in 0:max_days) {
    for (cand in c(d - k, d + k)) {
      p <- files[as.character(cand)]
      if (!is.na(p)) return(unname(p))
    }
  }
  stop("No file within ", max_days, " days of ", as.Date(d))
}

# Raw stored values (no scale/fill decoding, like mask_and_scale=False), as an [x, y]
# matrix in file order with the length-1 time dim dropped.
read_raw <- function(path, var) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  v <- ncvar_get(nc, var, raw_datavals = TRUE, collapse_degen = TRUE)
  dv <- nc$var[[var]]
  dn <- vapply(dv$dim, function(d) d$name, "")[dv$varsize > 1]
  if (identical(dn[1], "y")) v <- t(v)
  # Unsigned bytes read through a signed type come back negative (255 -> -1).
  neg <- !is.na(v) & v < 0
  if (any(neg)) v[neg] <- v[neg] + 256
  v
}


# ---------------------------------------------------------------------------
# Projection helpers
# ---------------------------------------------------------------------------
# Minimal CF grid_mapping -> PROJ for polar stereographic (stands in for pyproj's CRS.from_cf).
cf_polar_stereo <- function(a) {
  p <- c("+proj=stere",
         paste0("+lat_0=", a$latitude_of_projection_origin),
         paste0("+lon_0=", a$straight_vertical_longitude_from_pole %||% a$longitude_of_projection_origin %||% 0))
  if (!is.null(a$standard_parallel)) {
    p <- c(p, paste0("+lat_ts=", a$standard_parallel))
  } else if (!is.null(a$scale_factor_at_projection_origin)) {
    p <- c(p, paste0("+k=", a$scale_factor_at_projection_origin))
  }
  p <- c(p, paste0("+x_0=", a$false_easting %||% 0), paste0("+y_0=", a$false_northing %||% 0))
  if (!is.null(a$semi_major_axis)) {
    if (!is.null(a$inverse_flattening) && a$inverse_flattening > 0) {
      p <- c(p, paste0("+a=", a$semi_major_axis), paste0("+rf=", a$inverse_flattening))
    } else if (!is.null(a$semi_minor_axis)) {
      p <- c(p, paste0("+a=", a$semi_major_axis), paste0("+b=", a$semi_minor_axis))
    } else {
      p <- c(p, paste0("+R=", a$semi_major_axis))
    }
  } else if (!is.null(a$earth_radius)) {
    p <- c(p, paste0("+R=", a$earth_radius))
  } else {
    p <- c(p, "+datum=WGS84")
  }
  st_crs(paste(c(p, "+units=m", "+no_defs"), collapse = " "))
}

grid_crs <- function(nc, crs_var) {
  a <- ncatt_get(nc, crs_var)
  if (identical(a$grid_mapping_name, "polar_stereographic")) {
    crs <- tryCatch(cf_polar_stereo(a), error = function(e) NULL)
    if (!is.null(crs)) return(crs)
  }
  for (key in c("crs_wkt", "spatial_ref", "proj4text", "proj4_string", "proj4")) {
    if (!is.null(a[[key]])) return(st_crs(a[[key]]))
  }
  stop("Could not read a CRS from '", crs_var, "' attributes: ", paste(names(a), collapse = ", "))
}

# Fractional (col, row) position of each lat/lon inside the file's x/y grid (0-based, file order).
grid_fraction <- function(sample_path, crs_var, lat, lon) {
  nc <- nc_open(sample_path)
  on.exit(nc_close(nc))
  crs <- grid_crs(nc, crs_var)
  xs <- as.numeric(nc$dim$x$vals)
  ys <- as.numeric(nc$dim$y$vals)
  xy <- sf_project("EPSG:4326", crs$wkt, cbind(lon, lat), keep = TRUE, warn = FALSE)
  fx <- (xy[, 1] - xs[1]) / (xs[2] - xs[1])
  fy <- (xy[, 2] - ys[1]) / (ys[2] - ys[1])
  fx[!is.finite(fx)] <- -9
  fy[!is.finite(fy)] <- -9
  list(fx = fx, fy = fy, nx = length(xs), ny = length(ys))
}


# ---------------------------------------------------------------------------
# Globe geometry (computed once)
# ---------------------------------------------------------------------------
load_blue_marble <- function(month) {
  p <- sort(list.files(BM_DIR, pattern = sprintf("^world\\.2004%02d.*\\.jpg$", month), full.names = TRUE))
  p <- p[!startsWith(basename(p), "._")]
  if (length(p) == 0) stop(sprintf("No Blue Marble image for month %02d in %s", month, BM_DIR))
  readJPEG(p[1])   # [h, w, 3] in 0-1
}

# Anti-aliased lat circles (every 10 deg) and meridians (every 30 deg), in exact px distance.
graticule <- function(lat, rc, lon, hemi) {
  sgn <- if (hemi == "north") 1 else -1
  d_lat <- rep(Inf, length(rc))
  for (k in seq(10, 80, by = 10)) d_lat <- pmin(d_lat, R * abs(rc - cos(k * pi / 180)))
  dlon <- ((lon + 15) %% 30) - 15
  d_lon <- R * rc * abs(sin(dlon * pi / 180))
  d_lon[sgn * lat > 80] <- Inf                    # stop meridians before the pole
  pmin(pmax(1 - pmin(d_lat, d_lon) / 0.8, 0), 1)
}

build_globe <- function(hemi, bm_month, ims_sample = NULL, sic_sample = NULL) {
  cx <- CENTERS[[hemi]][1]; cy <- CENTERS[[hemi]][2]
  ys_seq <- (cy - R):(cy + R - 1)
  xs_seq <- (cx - R):(cx + R - 1)
  ys <- rep(ys_seq, times = length(xs_seq))
  xs <- rep(xs_seq, each = length(ys_seq))
  u <- (xs + 0.5 - cx) / R
  v <- (cy - (ys + 0.5)) / R
  r <- sqrt(u^2 + v^2)
  inside <- r < 1
  ys <- ys[inside]; xs <- xs[inside]; u <- u[inside]; v <- v[inside]; r <- r[inside]
  rc <- pmin(r, 1)
  
  if (hemi == "north") {           # Greenwich meridian pointing down, 90E to the right
    lat <- acos(rc) * 180 / pi
    lon <- atan2(u, -v) * 180 / pi
  } else {                         # Greenwich meridian pointing up, 90E to the right
    lat <- -acos(rc) * 180 / pi
    lon <- atan2(u, v) * 180 / pi
  }
  
  lin <- ys + xs * H + 1           # column-major index into an H x W plane (0-based px)
  g <- list(
    n = length(lat), lat = lat, lon = lon,
    lin3 = c(lin, lin + H * W, lin + 2 * H * W),   # R, G, B planes - matches n x 3 column order
    grid = graticule(lat, rc, lon, hemi) * GRID_ALPHA,
    edge = pmin(pmax((1 - r) * R, 0), 1),
    haze = r^HAZE_POWER * HAZE_STRENGTH
  )
  
  img <- load_blue_marble(bm_month)
  h <- dim(img)[1]; w <- dim(img)[2]
  col <- pmin(pmax(floor((lon + 180) / 360 * w), 0), w - 1)
  row <- pmin(pmax(floor((90 - lat) / 180 * h), 0), h - 1)
  g$base <- round(sapply(1:3, function(k) img[cbind(row + 1, col + 1, k)]) * 255)
  rm(img); gc(verbose = FALSE)
  
  if (hemi == "north") {
    gf <- grid_fraction(ims_sample, "projection", lat, lon)
    ix <- round(gf$fx); iy <- round(gf$fy)      # round-half-even, same as np.rint
    ok <- ix >= 0 & ix < gf$nx & iy >= 0 & iy < gf$ny
    g$ims_sel <- which(ok)
    g$ims_idx <- cbind(ix[ok] + 1, iy[ok] + 1)   # [x, y] matrix index
  } else {
    sel <- which(lat < -30)
    gf <- grid_fraction(sic_sample, "crs", lat[sel], lon[sel])
    x0 <- floor(gf$fx); y0 <- floor(gf$fy)
    ok <- x0 >= 0 & x0 < gf$nx - 1 & y0 >= 0 & y0 < gf$ny - 1
    g$sic_sel <- sel[ok]
    g$sic_x0 <- x0[ok] + 1; g$sic_y0 <- y0[ok] + 1   # 1-based
    g$sic_wx <- gf$fx[ok] - x0[ok]
    g$sic_wy <- gf$fy[ok] - y0[ok]
  }
  g
}


# ---------------------------------------------------------------------------
# Per-frame compositing
# ---------------------------------------------------------------------------
blend <- function(rgb, color, alpha) rgb * (1 - alpha) + outer(alpha, color)

composite <- function(g, canvas, ice_alpha = NULL, snow_alpha = NULL) {
  rgb <- g$base
  if (!is.null(snow_alpha)) rgb <- blend(rgb, SNOW_RGB, snow_alpha)
  if (!is.null(ice_alpha))  rgb <- blend(rgb, ICE_RGB, ice_alpha)
  rgb <- blend(rgb, GRID_RGB, g$grid)
  rgb <- blend(rgb, HAZE_RGB, g$haze)
  rgb <- rgb * g$edge
  canvas[g$lin3] <- floor(pmin(pmax(rgb, 0), 255)) / 255   # uint8 truncation, then 0-1 for rasterImage
  canvas
}

north_layers <- function(g, ims) {
  vals <- ims[g$ims_idx]
  snow <- ice <- numeric(g$n)
  snow[g$ims_sel[vals == 4]] <- SNOW_ALPHA
  ice[g$ims_sel[vals == 3]]  <- ICE_ALPHA
  list(ice = ice, snow = snow)
}

south_layer <- function(g, sic) {
  ok <- sic <= 100                                  # 255 = land / fill
  valid <- ok * 1
  A <- pmin(pmax((ifelse(ok, sic, 0) - 15) / (SIC_FULL - 15), 0), 1)^SIC_GAMMA
  wx <- g$sic_wx; wy <- g$sic_wy
  # Bilinear over ocean cells only, so land cells don't thin the ice along the coast
  num <- den <- numeric(length(wx))
  for (o in list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))) {    # (dx, dy)
    w <- (if (o[1]) wx else 1 - wx) * (if (o[2]) wy else 1 - wy)
    idx <- cbind(g$sic_x0 + o[1], g$sic_y0 + o[2])
    vw <- valid[idx] * w
    num <- num + A[idx] * vw
    den <- den + vw
  }
  ice <- numeric(g$n)
  ice[g$sic_sel] <- ifelse(den > 0, num / den, 0) * ICE_ALPHA
  ice
}


# ---------------------------------------------------------------------------
# Frame output / overlays
# ---------------------------------------------------------------------------
format_date <- function(d) sprintf("%s %d, %d", format(d, "%b"), as.integer(format(d, "%d")),
                                   as.integer(format(d, "%Y")))      # e.g. Aug 5, 2019

# res = 100 makes cex = pt / 12 match matplotlib font sizes at dpi = 100.
write_frame <- function(path, canvas, d) {
  png(path, width = W, height = H, res = 100, bg = "black", type = "cairo", family = FONT_FAMILY)
  on.exit(dev.off())
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(c(0, W), c(0, H))
  rasterImage(canvas, 0, 0, W, H, interpolate = FALSE)
  # Python layout uses top-down y; flip with H - y.
  text(W / 2, H - 58, "Snow Cover & Sea Ice", col = "white", cex = 30 / 12, font = 2)
  text(CENTERS$north[1], H - 112, "Northern Hemisphere", col = "white", cex = 20 / 12)
  text(CENTERS$south[1], H - 112, "Southern Hemisphere", col = "white", cex = 20 / 12)
  text(DATE_XY[1], H - DATE_XY[2], format_date(d), col = "white", cex = 24 / 12, adj = c(1, 0.5))
}

# Confirms the encoded video is real and complete before it leaves local disk.
.verify_video <- function(path, expected_frames) {
  if (!file.exists(path) || file.size(path) == 0) stop("Video missing or empty: ", path)
  info <- av_media_info(path)
  n <- info$video$frames
  message(sprintf("Verified %s: %s frames, %.1f s", basename(path), n, info$duration))
  if (is.null(n) || is.na(n) || n < expected_frames) {
    stop("Frame count lower than expected - not moving to drive.")
  }
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  if (!dir.exists(DRIVE)) stop(DRIVE, " is not mounted.")
  dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  ims_files <- index_files(IMS_DIR, "ims([0-9]{7})_", "%Y%j")
  sic_files <- index_files(SIC_DIR, "sic_pss25_([0-9]{8})_", "%Y%m%d")
  message(sprintf("Indexed %d IMS files, %d G02202 south files", length(ims_files), length(sic_files)))
  
  message("Precomputing globe geometry ...")
  north <- build_globe("north", NORTH_BM_MONTH, ims_sample = ims_files[[1]])
  south <- build_globe("south", SOUTH_BM_MONTH, sic_sample = sic_files[[1]])
  
  canvas <- array(0, dim = c(H, W, 3))
  
  draw <- function(d, canvas) {
    ims <- read_raw(nearest_file(ims_files, d), IMS_VAR)
    nl <- north_layers(north, ims)
    rm(ims)   # full 6144 x 6144 IMS grid; drop it before the next read
    canvas <- composite(north, canvas, ice_alpha = nl$ice, snow_alpha = nl$snow)
    ice_s <- south_layer(south, read_raw(nearest_file(sic_files, d), SIC_VAR))
    composite(south, canvas, ice_alpha = ice_s)
  }
  
  if (PREVIEW_ONLY) {
    out <- file.path(SCRIPT_DIR, sprintf("preview_%s_R.png", PREVIEW_DATE))
    write_frame(out, draw(PREVIEW_DATE, canvas), PREVIEW_DATE)
    message("Preview written: ", out)
    return(invisible())
  }
  
  dates <- seq(START, END, by = FRAME_STEP_DAYS)
  message(sprintf("Rendering %d frames at %d fps (~%.0f s)", length(dates), FPS, length(dates) / FPS))
  dir.create(FRAME_DIR, recursive = TRUE, showWarnings = FALSE)
  frames <- file.path(FRAME_DIR, sprintf("frame_%05d.png", seq_along(dates)))
  
  pb <- txtProgressBar(min = 0, max = length(dates), style = 3)
  for (i in seq_along(dates)) {
    canvas <- draw(dates[i], canvas)
    write_frame(frames[i], canvas, dates[i])
    setTxtProgressBar(pb, i)
    if (i %% 20 == 0) gc(verbose = FALSE)
  }
  close(pb)
  
  message("Encoding video...")
  av_encode_video(frames, output = LOCAL_OUT, framerate = FPS, codec = "libx264", verbose = FALSE)
  .verify_video(LOCAL_OUT, length(dates))
  
  # file.rename() fails across filesystems (APFS -> exFAT); copy, confirm, then delete.
  dir.create(dirname(FINAL_OUT), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(LOCAL_OUT, FINAL_OUT, overwrite = TRUE) ||
      file.size(FINAL_OUT) != file.size(LOCAL_OUT)) {
    stop("Copy to ", FINAL_OUT, " failed or is incomplete; local copy kept at ", LOCAL_OUT)
  }
  unlink(LOCAL_OUT)
  unlink(FRAME_DIR, recursive = TRUE)
  message("Moved to ", FINAL_OUT)
}

main()