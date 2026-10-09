# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations
# Library by Professor John Michael Wallace.
#
# Script: animate_fires_nicolosi.R
#
# Description: Generates the fire carbon emissions time-series animation from
#              GFED5.1 weekly totals (March 2003 - January 2022), rendered as two
#              Nicolosi globular hemispheres (west centered on 90W, east on 90E) in
#              the style of the NASA SVS "Fires - a global perspective" animation.
#
# Note: For the Plate Carree, Robinson, and Foucaut projections, see the other
#       scripts in the fires scripts folder.
#
# Run:
#   1. RENDER_MODE <- "test"  -> writes a few PNG frames for review (fast)
#   2. RENDER_MODE <- "full"  -> caffeinate -i Rscript animate_fires_nicolosi.R.


Sys.setenv(PROJ_DEBUG = "0")   # silence PROJ notes
suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(png)
  library(av)
})

# ---------------------------------------------------------------- settings
RENDER_MODE <- "full"            # "test" or "full"
TEST_DATES  <- as.Date(c("2004-03-06", "2010-08-14", "2015-09-26", "2019-12-28"))

WORK        <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "fires")   # local APFS
WEEKLY      <- "/Volumes/CLIMATEDATA/fires/GFED5.1_weekly"
BASEMAP_DIR <- "/Volumes/CLIMATEDATA/fires/basemap"
STAGING     <- file.path(WORK, "_staging")
TEST_DIR    <- file.path(WORK, "test_frames")
FRAME_DIR   <- "/Volumes/CLIMATEDATA/_frames_tmp_fires_nicolosi_R"
BG_CACHE    <- file.path(WORK, "fires_background_nicolosi_R.png")
FINAL_VIDEO <- "/Volumes/CLIMATEDATA/fires_2003_2022_nicolosi_R.mp4"

WIDTH <- 1920; HEIGHT <- 1080; RES <- 100   # px; RES makes font sizes match the Python version
FPS   <- 12                                 # 987 weeks -> ~82 s (NASA original runs 1:24)
EARTH_RADIUS <- 6371007.2                   # m (sphere)
HEMI_CENTERS <- c(-90, 90)                  # central longitudes: western, eastern hemisphere
HEMI_GAP     <- 24                          # px between the two circles
MAP_MARGIN_X <- 24                          # min px on each side
MAP_MARGIN_Y <- 24                          # min px top and bottom
SUPERSAMPLE  <- 4                           # lon/lat samples per data cell, per axis, for the pixel table
LANDMASK_RES <- 0.05                        # degrees; resolution of the rasterized land mask

# Basemap colors (NASA SVS look)
FRAME_COLOR   <- "#1c1e22"                  # outside the map outline
OCEAN_COLOR   <- "#4d5868"
LAND_DARK     <- 0.10; LAND_LIGHT <- 0.36   # gray range for shaded relief (0 = black, 1 = white)
LAND_FLAT     <- "#2e2e2e"                  # used only if the relief raster is unavailable
BORDER_COLOR  <- adjustcolor("#d9d9d9", alpha.f = 0.6)
OUTLINE_COLOR <- "#8b919b"
RELIEF_URLS   <- c(
  "https://naciscdn.org/naturalearth/50m/raster/GRAY_50M_SR.zip",
  "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/GRAY_50M_SR.zip"
)

# Data colors: green -> orange -> white, 0-20 g C m-2 week-1, saturating above 20
VMIN <- 0; VMAX <- 20
CMAP_NODES <- c("#5e8c3a", "#a3b257", "#e5b260", "#f0d5a4", "#f3f1ec")   # at 0, .25, .5, .75, 1
ALPHA_MIN  <- 0.05; ALPHA_FULL <- 1.0      # transparent below ALPHA_MIN, opaque above ALPHA_FULL
TEXT_COLOR <- "#e6e6e6"
HALO_COLOR <- "#1c1e22"                     # outline behind legend text, for legibility over land
CREDIT     <- "Data: GFED5.1 (van der Werf et al.)  |  SDSU Climate Informatics Lab"

LUT_STEP <- 0.02                            # g C m-2 week-1 per color-table bin

invisible(Sys.setlocale("LC_TIME", "C"))    # English month abbreviations
PNG_TYPE <- if (capabilities("cairo")) "cairo" else NULL


# ---------------------------------------------------------------- projection geometry
hemi_crs <- function(lon0) sprintf("+proj=nicol +lon_0=%g +R=%.1f +units=m +no_defs", lon0, EARTH_RADIUS)

map_geometry <- function() {
  # Two circles side by side, as large as fit inside the margins, centered
  r_m <- sf_project("EPSG:4326", hemi_crs(HEMI_CENTERS[2]), cbind(HEMI_CENTERS[2] + 90, 0))[1, 1]
  d <- min((WIDTH - 2 * MAP_MARGIN_X - HEMI_GAP) / 2, HEIGHT - 2 * MAP_MARGIN_Y)
  x_start <- (WIDTH - (2 * d + HEMI_GAP)) / 2
  list(r_px = d / 2, scale = (d / 2) / r_m,
       cx = c(x_start + d / 2, x_start + 1.5 * d + HEMI_GAP), cy = c(HEIGHT / 2, HEIGHT / 2))
}

hemisphere <- function(lon) ifelse(lon >= 0, 2L, 1L)   # 1 = western, 2 = eastern

lonlat_to_px <- function(G, lon, lat) {
  # Continuous pixel coordinates (col, row from the top-left corner)
  col <- rep(NA_real_, length(lon)); row <- col
  h <- hemisphere(lon)
  for (k in 1:2) {
    m <- which(h == k)
    if (length(m)) {
      xy <- sf_project("EPSG:4326", hemi_crs(HEMI_CENTERS[k]), cbind(lon[m], lat[m]),
                       keep = TRUE, warn = FALSE)
      col[m] <- G$cx[k] + xy[, 1] * G$scale
      row[m] <- G$cy[k] - xy[, 2] * G$scale
    }
  }
  list(col = col, row = row)
}

inside_mask <- function(G) {
  cc <- seq_len(WIDTH) - 0.5; rr <- seq_len(HEIGHT) - 0.5
  inside <- matrix(FALSE, HEIGHT, WIDTH)
  for (k in 1:2) inside <- inside | outer((rr - G$cy[k])^2, (cc - G$cx[k])^2, "+") <= G$r_px^2
  inside
}

shift <- function(m, dr, dc) {
  # Shift a matrix by (dr, dc) cells, filling with NA (no wrap-around)
  out <- matrix(NA_real_, nrow(m), ncol(m))
  rs <- max(1, 1 - dr):min(nrow(m), nrow(m) - dr)
  cs <- max(1, 1 - dc):min(ncol(m), ncol(m) - dc)
  out[rs + dr, cs + dc] <- m[rs, cs]
  out
}

pixel_lonlat <- function(G, cell_deg) {
  # lon/lat for every pixel inside the circles, by forward-projecting a dense sample grid.
  # Returns HEIGHT x WIDTH matrices (NA outside the map).
  step <- cell_deg / SUPERSAMPLE
  lons <- -180 + (seq_len(round(360 / step)) - 0.5) * step
  lats <- -90 + (seq_len(round(180 / step)) - 0.5) * step
  pix_lon <- rep(NA_real_, WIDTH * HEIGHT); pix_lat <- pix_lon

  bands <- split(seq_along(lats), ceiling(seq_along(lats) / 160))
  pb <- txtProgressBar(min = 0, max = length(bands), style = 3)
  for (b in seq_along(bands)) {
    LO <- rep(lons, times = length(bands[[b]]))
    LA <- rep(lats[bands[[b]]], each = length(lons))
    p  <- lonlat_to_px(G, LO, LA)
    ok <- is.finite(p$col) & is.finite(p$row)
    c0 <- floor(p$col[ok]); r0 <- floor(p$row[ok])
    inb <- c0 >= 0 & c0 < WIDTH & r0 >= 0 & r0 < HEIGHT
    idx <- r0[inb] * WIDTH + c0[inb] + 1                # row-major, 1-based
    pix_lon[idx] <- LO[ok][inb]
    pix_lat[idx] <- LA[ok][inb]
    setTxtProgressBar(pb, b)
  }
  close(pb)
  pix_lon <- matrix(pix_lon, nrow = HEIGHT, byrow = TRUE)
  pix_lat <- matrix(pix_lat, nrow = HEIGHT, byrow = TRUE)

  # Fill the few pixels inside the circles that no sample landed in
  inside <- inside_mask(G)
  holes0 <- sum(inside & is.na(pix_lon))
  for (iter in 1:20) {
    holes <- inside & is.na(pix_lon)
    if (!any(holes)) break
    for (s in list(c(0, 1), c(0, -1), c(1, 0), c(-1, 0))) {
      src_lon <- shift(pix_lon, s[1], s[2]); src_lat <- shift(pix_lat, s[1], s[2])
      take <- holes & !is.na(src_lon)
      pix_lon[take] <- src_lon[take]; pix_lat[take] <- src_lat[take]
      holes <- holes & !take
    }
  }
  message(sprintf("Pixel table: %s map pixels (%s gap pixels filled, %d left)",
                  format(sum(!is.na(pix_lon)), big.mark = ","), format(holes0, big.mark = ","),
                  sum(inside & is.na(pix_lon))))
  list(lon = pix_lon, lat = pix_lat)
}


# ---------------------------------------------------------------- helpers
move_file <- function(from, to) {
  # Cross-filesystem move (local APFS -> exFAT): copy, check size, then delete
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok || file.size(to) != file.size(from)) stop("Copy to ", to, " failed; left at ", from)
  unlink(from)
}

open_png <- function(path, width = WIDTH, height = HEIGHT, bg = FRAME_COLOR,
                     xlim = c(0, WIDTH), ylim = c(0, HEIGHT)) {
  args <- list(filename = path, width = width, height = height, res = RES, bg = bg)
  if (!is.null(PNG_TYPE)) args$type <- PNG_TYPE
  do.call(png, args)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = xlim, ylim = ylim)       # default: 1 user unit = 1 px, y up
}

nc_dates <- function(nc) {
  vals  <- ncvar_get(nc, "time")
  units <- ncatt_get(nc, "time", "units")$value          # e.g. "days since 2003-03-01"
  step  <- tolower(strsplit(units, " ")[[1]][1])
  origin <- as.POSIXct(sub(".*since ", "", units), tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"))
  secs <- switch(step, days = 86400, hours = 3600, minutes = 60, seconds = 1,
                 stop("Unsupported time units: ", units))
  as.Date(origin + vals * secs)
}


# ---------------------------------------------------------------- basemap
load_relief <- function() {
  # Natural Earth 1:50m gray shaded relief (shared with the Python scripts) as a
  # north-up matrix (5400 x 2700 after 2x averaging); NULL if unavailable
  dir.create(BASEMAP_DIR, recursive = TRUE, showWarnings = FALSE)
  tif <- file.path(BASEMAP_DIR, "GRAY_50M_SR.tif")
  if (!file.exists(tif)) {
    for (url in RELIEF_URLS) {
      message("Downloading shaded relief: ", url)
      zf <- tempfile(fileext = ".zip")
      ok <- tryCatch({ download.file(url, zf, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) { message("  failed: ", conditionMessage(e)); FALSE })
      if (ok) {
        members <- unzip(zf, list = TRUE)$Name
        member  <- members[grepl("\\.tif$", members, ignore.case = TRUE) &
                             !grepl("(^|/)\\._|__MACOSX", members)][1]
        unzip(zf, files = member, exdir = BASEMAP_DIR, junkpaths = TRUE)
        file.rename(file.path(BASEMAP_DIR, basename(member)), tif)
        break
      }
    }
  }
  if (!file.exists(tif)) { message("Shaded relief unavailable; using flat land color."); return(NULL) }
  r <- suppressWarnings(rast(tif))[[1]]                  # the .tif has no embedded georeference
  ext(r) <- ext(-180, 180, -90, 90)
  as.matrix(aggregate(r, fact = 2, fun = "mean"), wide = TRUE)
}

ne_layer <- function(type, category) {
  # Natural Earth 1:50m vector layer, downloaded once into BASEMAP_DIR
  dir.create(BASEMAP_DIR, recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    ne_load(scale = 50, type = type, category = category, destdir = BASEMAP_DIR, returnclass = "sf"),
    error = function(e) ne_download(scale = 50, type = type, category = category,
                                    destdir = BASEMAP_DIR, load = TRUE, returnclass = "sf")
  )
}

land_fraction_raster <- function() {
  # Antialiased land mask (1 = land, 0 = ocean/lake) on a regular lon/lat grid, north up,
  # drawn with R's PNG device (the same idea as the matplotlib version)
  nx <- round(360 / LANDMASK_RES); ny <- round(180 / LANDMASK_RES)
  tmp <- tempfile(fileext = ".png")
  open_png(tmp, width = nx, height = ny, bg = "black", xlim = c(-180, 180), ylim = c(-90, 90))
  plot(st_geometry(ne_layer("land", "physical")), col = "white", border = NA, add = TRUE)
  plot(st_geometry(ne_layer("lakes", "physical")), col = "black", border = NA, add = TRUE)
  invisible(dev.off())
  mask <- readPNG(tmp)[, , 1]
  unlink(tmp)
  mask
}

sample_bilinear <- function(grid, lon, lat) {
  # Bilinear sample of a north-up global lon/lat grid at the given points
  ny <- nrow(grid); nx <- ncol(grid)
  fx <- (lon + 180) / 360 * nx - 0.5
  fy <- (90 - lat) / 180 * ny - 0.5
  x0 <- floor(fx); y0 <- floor(fy)
  tx <- fx - x0; ty <- fy - y0
  x0c <- (x0 %% nx) + 1; x1c <- ((x0 + 1) %% nx) + 1            # wrap in longitude
  y0c <- pmin(pmax(y0, 0), ny - 1) + 1; y1c <- pmin(pmax(y0 + 1, 0), ny - 1) + 1
  top <- grid[cbind(y0c, x0c)] * (1 - tx) + grid[cbind(y0c, x1c)] * tx
  bot <- grid[cbind(y1c, x0c)] * (1 - tx) + grid[cbind(y1c, x1c)] * tx
  top * (1 - ty) + bot * ty
}

border_polyline <- function(G) {
  # Country borders forward-projected to pixel coordinates (y up), with NA breaks
  # between lines and wherever a line crosses from one hemisphere to the other
  g <- st_geometry(ne_layer("admin_0_boundary_lines_land", "cultural"))
  g <- st_segmentize(st_set_crs(g, NA), dfMaxLength = 0.25)   # planar, in degrees
  g <- st_cast(g, "MULTILINESTRING")                           # uniform type -> columns X, Y, L1, L2
  xy <- st_coordinates(g)
  # One id per line part (st_coordinates gives L1, or L1 + L2 for multi-part lines)
  ids <- xy[, setdiff(colnames(xy), c("X", "Y")), drop = FALSE]
  grp <- do.call(paste, c(lapply(seq_len(ncol(ids)), function(j) ids[, j]), sep = "_"))
  stopifnot(length(grp) == nrow(xy))
  p <- lonlat_to_px(G, xy[, "X"], xy[, "Y"])
  brk <- c(TRUE, grp[-1] != grp[-length(grp)] |
             diff(hemisphere(xy[, "X"])) != 0)              # start of a new piece
  x <- p$col; y <- HEIGHT - p$row
  n <- length(x)
  out_x <- rep(NA_real_, n + sum(brk)); out_y <- out_x
  pos <- seq_len(n) + cumsum(brk)                          # leave an NA slot before each piece
  out_x[pos] <- x; out_y[pos] <- y
  list(x = out_x, y = out_y)
}

build_background <- function(G, px) {
  # Render the static basemap once; returns a (HEIGHT, WIDTH, 3) array in 0-1
  if (file.exists(BG_CACHE)) return(readPNG(BG_CACHE)[, , 1:3])

  message("Building background ...")
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  valid <- which(!is.na(t(px$lon)))                         # row-major pixel order
  lon <- t(px$lon)[valid]; lat <- t(px$lat)[valid]

  land <- sample_bilinear(land_fraction_raster(), lon, lat)
  relief <- load_relief()
  if (!is.null(relief)) {
    ny <- nrow(relief); nx <- ncol(relief)
    ir <- pmin(pmax(floor((90 - lat) / 180 * ny) + 1, 1), ny)
    ic <- pmin(pmax(floor((lon + 180) / 360 * nx) + 1, 1), nx)
    v <- relief[cbind(ir, ic)]; v[is.na(v)] <- 0
    land_rgb <- matrix(rep(LAND_DARK + (LAND_LIGHT - LAND_DARK) * v / 255, 3), ncol = 3)
  } else {
    land_rgb <- matrix(col2rgb(LAND_FLAT) / 255, nrow = length(lon), ncol = 3, byrow = TRUE)
  }
  ocean_rgb <- matrix(col2rgb(OCEAN_COLOR) / 255, nrow = length(lon), ncol = 3, byrow = TRUE)
  mix <- ocean_rgb * (1 - land) + land_rgb * land
  img <- rep(FRAME_COLOR, WIDTH * HEIGHT)
  img[valid] <- rgb(mix[, 1], mix[, 2], mix[, 3])

  borders <- border_polyline(G)
  open_png(BG_CACHE)
  rasterImage(as.raster(matrix(img, nrow = HEIGHT, byrow = TRUE)), 0, 0, WIDTH, HEIGHT,
              interpolate = FALSE)
  lines(borders$x, borders$y, col = BORDER_COLOR, lwd = 0.6)
  a <- seq(0, 2 * pi, length.out = 721)
  for (k in 1:2) {
    polygon(G$cx[k] + G$r_px * cos(a), (HEIGHT - G$cy[k]) + G$r_px * sin(a),
            border = OUTLINE_COLOR, col = NA, lwd = 1.05)
  }
  invisible(dev.off())

  bg <- readPNG(BG_CACHE)[, , 1:3]
  if (any(dim(bg)[1:2] != c(HEIGHT, WIDTH))) stop("Background has wrong size: ", paste(dim(bg), collapse = "x"))
  bg
}


# ---------------------------------------------------------------- colors and frames
make_lut <- function() {
  # One color (with alpha) per LUT_STEP bin from 0 to VMAX; the last bin is the >= VMAX color
  v <- seq(0, VMAX, by = LUT_STEP)
  rgb_v <- colorRamp(CMAP_NODES)(pmin(pmax((v - VMIN) / (VMAX - VMIN), 0), 1))
  alpha <- pmin(pmax((v - ALPHA_MIN) / (ALPHA_FULL - ALPHA_MIN), 0), 1)
  rgb(rgb_v[, 1], rgb_v[, 2], rgb_v[, 3], alpha * 255, maxColorValue = 255)
}

build_lookup <- function(px, lat, lon) {
  # Map pixels (row-major) and their source cell as a column-major index into
  # ncdf4's [lon, lat] matrix
  plon <- as.vector(t(px$lon)); plat <- as.vector(t(px$lat))
  pix <- which(!is.na(plon))
  dlat <- lat[2] - lat[1]; dlon <- lon[2] - lon[1]
  ilat <- pmin(pmax(round((plat[pix] - lat[1]) / dlat) + 1, 1), length(lat))
  ilon <- pmin(pmax(round((plon[pix] - lon[1]) / dlon) + 1, 1), length(lon))
  list(pix = pix, src = ilon + (ilat - 1) * length(lon))
}

fire_raster <- function(field, lut, lk) {
  # field: [lon, lat] matrix from ncdf4 -> full-frame raster (transparent outside the map)
  idx <- pmin(floor(as.vector(field) / LUT_STEP) + 1, length(lut))
  idx[!is.finite(idx) | idx < 1] <- 1
  buf <- rep("#00000000", WIDTH * HEIGHT)
  buf[lk$pix] <- lut[idx[lk$src]]
  as.raster(matrix(buf, nrow = HEIGHT, byrow = TRUE))
}

halo_text <- function(x, y, labels, ..., r = 1.5) {
  # Text with a dark outline (R's equivalent of the Python path-effect halo); 1 user unit = 1 px
  for (a in seq(0, 2 * pi, length.out = 9)[-9]) {
    text(x + r * cos(a), y + r * sin(a), labels, col = HALO_COLOR, ...)
  }
  text(x, y, labels, col = TEXT_COLOR, ...)
}

draw_legend <- function(date_label) {
  # Bottom center, in the open space between the two hemispheres (figure fractions as in Python)
  w <- 0.20; x0 <- 0.5 - w / 2; y0 <- 0.075; h <- 0.022
  X <- function(f) f * WIDTH; Y <- function(f) f * HEIGHT
  bar_w <- w / 1.05                                  # leave room for the "extend" triangle
  n <- 256
  xs <- seq(x0, x0 + bar_w, length.out = n + 1)
  cols <- colorRampPalette(CMAP_NODES)(n)
  rect(X(xs[-(n + 1)]), Y(y0), X(xs[-1]), Y(y0 + h), col = cols, border = NA)
  polygon(X(c(x0 + bar_w, x0 + w, x0 + bar_w)), Y(c(y0, y0 + h / 2, y0 + h)),
          col = tail(CMAP_NODES, 1), border = NA)

  ticks <- seq(0, 20, by = 5)
  tx <- x0 + bar_w * (ticks - VMIN) / (VMAX - VMIN)
  segments(X(tx), Y(y0), X(tx), Y(y0 - 0.006), col = TEXT_COLOR, lwd = 1)
  halo_text(X(tx), Y(y0 - 0.010), ticks, adj = c(0.5, 1), cex = 14 / 12)
  halo_text(X(0.5), Y(y0 - 0.044), "grams of carbon per square meter per week",
            adj = c(0.5, 1), cex = 14 / 12)

  halo_text(X(0.5), Y(0.150), date_label, adj = c(0.5, 0), cex = 22 / 12, font = 2)
  halo_text(X(0.5), Y(0.112), "Carbon Emissions", adj = c(0.5, 0), cex = 18 / 12)
  text(X(0.99), Y(0.012), CREDIT, adj = c(1, 0), col = "#bdbdbd", cex = 11 / 12)
}

render_frame <- function(path, bg_raster, fire, date) {
  open_png(path)
  rasterImage(bg_raster, 0, 0, WIDTH, HEIGHT, interpolate = FALSE)
  rasterImage(fire, 0, 0, WIDTH, HEIGHT, interpolate = FALSE)
  draw_legend(format(date, "%b  %Y"))
  invisible(dev.off())
}


# ---------------------------------------------------------------- data index
week_index <- function() {
  files <- list.files(WEEKLY, pattern = "^GFED5\\.1_C_weekly_\\d{4}\\.nc$", full.names = TRUE)
  if (length(files) == 0) stop("No weekly files in ", WEEKLY)
  do.call(rbind, lapply(sort(files), function(f) {
    nc <- nc_open(f); on.exit(nc_close(nc))
    d <- nc_dates(nc)
    data.frame(file = f, k = seq_along(d), date = d, stringsAsFactors = FALSE)
  }))
}

read_week <- function(nc, k) {
  ncvar_get(nc, "C_weekly", start = c(1, 1, k), count = c(-1, -1, 1))   # [lon, lat]
}


# ---------------------------------------------------------------- main
main <- function() {
  if (!dir.exists(WEEKLY)) stop("Weekly data folder not found: ", WEEKLY)
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  dir.create(STAGING, recursive = TRUE, showWarnings = FALSE)

  index <- week_index()
  message(nrow(index), " weekly frames: ", index$date[1], " to ", index$date[nrow(index)])
  nc0 <- nc_open(index$file[1])
  lat <- ncvar_get(nc0, "lat"); lon <- ncvar_get(nc0, "lon")
  nc_close(nc0)

  G  <- map_geometry()
  px <- pixel_lonlat(G, abs(lon[2] - lon[1]))
  bg <- build_background(G, px)
  bg_raster <- as.raster(bg)
  lk  <- build_lookup(px, lat, lon)
  rm(px); invisible(gc())
  lut <- make_lut()

  if (RENDER_MODE == "test") {
    dir.create(TEST_DIR, recursive = TRUE, showWarnings = FALSE)
    for (d in as.list(TEST_DATES)) {
      i <- which.min(abs(as.numeric(index$date - d)))
      nc <- nc_open(index$file[i])
      field <- read_week(nc, index$k[i])
      nc_close(nc)
      out <- file.path(TEST_DIR, sprintf("fires_nicolosi_R_%s.png", index$date[i]))
      render_frame(out, bg_raster, fire_raster(field, lut, lk), index$date[i])
      message("  wrote ", out)
    }
    message("Review the PNGs, then set RENDER_MODE <- \"full\".")
    return(invisible())
  }

  # Full render: frames -> FRAME_DIR (resumable), then encode
  dir.create(FRAME_DIR, recursive = TRUE, showWarnings = FALSE)
  n <- nrow(index)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  current_file <- ""; nc <- NULL
  for (i in seq_len(n)) {
    out <- file.path(FRAME_DIR, sprintf("frame_%04d.png", i))
    if (!file.exists(out) || file.size(out) == 0) {
      if (index$file[i] != current_file) {
        if (!is.null(nc)) nc_close(nc)
        nc <- nc_open(index$file[i]); current_file <- index$file[i]
      }
      render_frame(out, bg_raster, fire_raster(read_week(nc, index$k[i]), lut, lk), index$date[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  if (!is.null(nc)) nc_close(nc)

  frames <- file.path(FRAME_DIR, sprintf("frame_%04d.png", seq_len(n)))   # exact list: no ._ files
  tmp <- file.path(STAGING, basename(FINAL_VIDEO))
  message("Encoding ", n, " frames ...")
  av_encode_video(frames, output = tmp, framerate = FPS, codec = "libx264",
                  vfilter = "format=yuv420p", verbose = FALSE)

  # av_media_info()$video$frames overestimates by one; duration x fps is exact
  nframes <- round(av_media_info(tmp)$duration * FPS)
  if (length(nframes) == 0 || nframes != n) stop("Video has ", nframes, " frames, expected ", n, "; left at ", tmp)
  move_file(tmp, FINAL_VIDEO)
  unlink(FRAME_DIR, recursive = TRUE)
  message("Done: ", FINAL_VIDEO, "  (", nframes, " frames, ", round(nframes / FPS, 1), " s)")
}

main()