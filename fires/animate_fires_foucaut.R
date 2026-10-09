# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations
# Library by Professor John Michael Wallace.
#
# Script: animate_fires_foucaut.R
#
# Description: Generates the fire carbon emissions time-series animation from
#              GFED5.1 weekly totals (March 2003 - January 2022), rendered in the
#              Foucaut projection in the style of the NASA SVS "Fires - a global
#              perspective" animation.
#
# Note: For the Plate Carree, Robinson, and Nicolosi projections, see the other
#       scripts in the fires scripts folder.
#
# Run:
#   1. RENDER_MODE <- "test"  -> writes a few PNG frames for review (fast)
#   2. RENDER_MODE <- "full"  -> caffeinate -i Rscript animate_fires_foucaut.R
# Frames are written to FRAME_DIR and skipped if already present, so an
# interrupted full render resumes where it stopped.
# If you change any basemap colors or the layout, delete BG_CACHE.
#
# Method: uses PROJ's plain `fouc` (it has a working inverse; `fouc_s` is not used)
# on the same sphere as the Python version. Every output pixel inside the Foucaut
# outline is inverse-projected to
# lon/lat once (with a forward round-trip check), giving a pixel -> grid-cell
# lookup. The shaded relief is sampled through the same lookup, and each frame
# is a single vector index into the weekly field, which is far faster than
# re-projecting 987 rasters.
#
# Input data: GFED5.1_C_weekly_YYYY.nc (variable C_weekly, g C m-2 week-1),
# available from the fire_climatology_v1.0 GitHub Release.

Sys.setenv(PROJ_DEBUG = "0")   # silence PROJ's "Invalid latitude" notes for pixels outside the outline
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
FRAME_DIR   <- "/Volumes/CLIMATEDATA/_frames_tmp_fires_foucaut_R"
BG_CACHE    <- file.path(WORK, "fires_background_foucaut_R.png")
FINAL_VIDEO <- "/Volumes/CLIMATEDATA/fires_2003_2022_foucaut_R.mp4"

WIDTH <- 1920; HEIGHT <- 1080; RES <- 100   # px; RES makes font sizes match the Python version
FPS   <- 12                                 # 987 weeks -> ~82 s (NASA original runs 1:24)
PROJ_CRS     <- "+proj=fouc +lon_0=0 +R=6378137 +units=m +no_defs"   # sphere, as in the Python version
MAP_MARGIN_X <- 24                          # min px on each side
MAP_MARGIN_Y <- 24                          # min px top and bottom

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


# ---------------------------------------------------------------- map layout
map_layout <- function() {
  # Equal-scale Foucaut map, as large as fits inside the margins, centered.
  # Returns the projected coordinates of the whole device and meters per pixel.
  ol <- outline_xy()
  xmin <- min(ol[, 1]); xmax <- max(ol[, 1])
  ymin <- min(ol[, 2]); ymax <- max(ol[, 2])
  mpp  <- max((xmax - xmin) / (WIDTH - 2 * MAP_MARGIN_X),
              (ymax - ymin) / (HEIGHT - 2 * MAP_MARGIN_Y))
  left_px <- (WIDTH - (xmax - xmin) / mpp) / 2
  top_px  <- (HEIGHT - (ymax - ymin) / mpp) / 2
  x_left <- xmin - left_px * mpp
  y_top  <- ymax + top_px * mpp
  list(mpp = mpp, x_left = x_left, x_right = x_left + WIDTH * mpp,
       y_top = y_top, y_bottom = y_top - HEIGHT * mpp)
}

pixel_lonlat <- function(L) {
  # lon/lat at the center of every pixel inside the map (row-major pixel order)
  X <- L$x_left + (seq_len(WIDTH) - 0.5) * L$mpp
  Y <- L$y_top - (seq_len(HEIGHT) - 0.5) * L$mpp
  pts <- cbind(rep(X, times = HEIGHT), rep(Y, each = WIDTH))
  ll <- suppressWarnings(sf_project(PROJ_CRS, "EPSG:4326", pts, keep = TRUE, warn = FALSE))
  fin <- which(is.finite(ll[, 1]) & is.finite(ll[, 2]))
  back <- suppressWarnings(sf_project("EPSG:4326", PROJ_CRS, ll[fin, , drop = FALSE],
                                      keep = TRUE, warn = FALSE))
  ok <- abs(back[, 1] - pts[fin, 1]) < 1000 & abs(back[, 2] - pts[fin, 2]) < 1000   # round trip
  ok[is.na(ok)] <- FALSE
  pix <- fin[ok]
  message(sprintf("Lookup: %s map pixels", format(length(pix), big.mark = ",")))
  list(pix = pix, lon = ll[pix, 1], lat = ll[pix, 2])
}

outline_xy <- function() {
  # The +/-180 meridians traced pole to pole (a hair inside, so PROJ keeps their sides)
  lat <- seq(-90, 90, length.out = 721)
  ll <- rbind(cbind(-180 + 1e-7, lat), cbind(180 - 1e-7, rev(lat)))
  xy <- sf_project("EPSG:4326", PROJ_CRS, ll)
  if (!all(is.finite(xy))) stop("Foucaut outline could not be traced (non-finite points).")
  xy
}


# ---------------------------------------------------------------- helpers
move_file <- function(from, to) {
  # Cross-filesystem move (local APFS -> exFAT): copy, check size, then delete
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok || file.size(to) != file.size(from)) stop("Copy to ", to, " failed; left at ", from)
  unlink(from)
}

open_png <- function(path, L) {
  args <- list(filename = path, width = WIDTH, height = HEIGHT, res = RES, bg = FRAME_COLOR)
  if (!is.null(PNG_TYPE)) args$type <- PNG_TYPE
  do.call(png, args)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(L$x_left, L$x_right), ylim = c(L$y_bottom, L$y_top))
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

to_map <- function(x) {
  # Densify in lon/lat (so straight edges along the dateline follow the curved
  # outline), then project to Foucaut
  g <- st_geometry(x)
  crs0 <- st_crs(g)
  g <- st_segmentize(st_set_crs(g, NA), dfMaxLength = 0.5)   # planar, in degrees (no lwgeom needed)
  st_transform(st_set_crs(g, crs0), PROJ_CRS)
}

build_background <- function(L, px) {
  # Render the static basemap once; returns a (HEIGHT, WIDTH, 3) array in 0-1
  if (file.exists(BG_CACHE)) return(readPNG(BG_CACHE)[, , 1:3])
  
  message("Building background ...")
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  ocean   <- to_map(ne_layer("ocean", "physical"))
  lakes   <- to_map(ne_layer("lakes", "physical"))
  borders <- to_map(ne_layer("admin_0_boundary_lines_land", "cultural"))
  relief  <- load_relief()
  
  # Raster layer in pixel space: shaded relief inside the outline, frame color outside
  img <- rep(FRAME_COLOR, WIDTH * HEIGHT)
  if (!is.null(relief)) {
    ny <- nrow(relief); nx <- ncol(relief)
    ir <- pmin(pmax(floor((90 - px$lat) / 180 * ny) + 1, 1), ny)
    ic <- pmin(pmax(floor((px$lon + 180) / 360 * nx) + 1, 1), nx)
    v  <- relief[cbind(ir, ic)]
    v[is.na(v)] <- 0
    img[px$pix] <- gray(LAND_DARK + (LAND_LIGHT - LAND_DARK) * v / 255)
  } else {
    img[px$pix] <- OCEAN_COLOR
  }
  
  open_png(BG_CACHE, L)
  rasterImage(as.raster(matrix(img, nrow = HEIGHT, byrow = TRUE)),
              L$x_left, L$y_bottom, L$x_right, L$y_top, interpolate = FALSE)
  if (is.null(relief)) plot(to_map(ne_layer("land", "physical")), col = LAND_FLAT, border = NA, add = TRUE)
  plot(ocean, col = OCEAN_COLOR, border = NA, add = TRUE)
  plot(lakes, col = OCEAN_COLOR, border = NA, add = TRUE)
  plot(borders, col = BORDER_COLOR, lwd = 0.6, add = TRUE)
  polygon(outline_xy(), border = OUTLINE_COLOR, col = NA, lwd = 1.05)
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
  # Source cell for each map pixel, as a column-major index into ncdf4's [lon, lat] matrix
  dlat <- lat[2] - lat[1]; dlon <- lon[2] - lon[1]
  ilat <- pmin(pmax(round((px$lat - lat[1]) / dlat) + 1, 1), length(lat))
  ilon <- pmin(pmax(round((px$lon - lon[1]) / dlon) + 1, 1), length(lon))
  list(pix = px$pix, src = ilon + (ilat - 1) * length(lon))
}

fire_raster <- function(field, lut, lk) {
  # field: [lon, lat] matrix from ncdf4 -> full-frame raster (transparent outside the map)
  idx <- pmin(floor(as.vector(field) / LUT_STEP) + 1, length(lut))
  idx[!is.finite(idx) | idx < 1] <- 1
  buf <- rep("#00000000", WIDTH * HEIGHT)
  buf[lk$pix] <- lut[idx[lk$src]]
  as.raster(matrix(buf, nrow = HEIGHT, byrow = TRUE))
}

ndc_x <- function(f) grconvertX(f, "ndc", "user")
ndc_y <- function(f) grconvertY(f, "ndc", "user")

halo_text <- function(x, y, labels, upp, ..., r = 1.5) {
  # Text with a dark outline (R's equivalent of the Python path-effect halo)
  for (a in seq(0, 2 * pi, length.out = 9)[-9]) {
    text(x + r * upp * cos(a), y + r * upp * sin(a), labels, col = HALO_COLOR, ...)
  }
  text(x, y, labels, col = TEXT_COLOR, ...)
}

draw_legend <- function(date_label, upp) {
  x0 <- 0.03; w <- 0.22; y0 <- 0.185; h <- 0.022
  bar_w <- w / 1.05                                  # leave room for the "extend" triangle
  n <- 256
  xs <- seq(x0, x0 + bar_w, length.out = n + 1)
  cols <- colorRampPalette(CMAP_NODES)(n)
  rect(ndc_x(xs[-(n + 1)]), ndc_y(y0), ndc_x(xs[-1]), ndc_y(y0 + h), col = cols, border = NA)
  polygon(ndc_x(c(x0 + bar_w, x0 + w, x0 + bar_w)), ndc_y(c(y0, y0 + h / 2, y0 + h)),
          col = tail(CMAP_NODES, 1), border = NA)
  
  ticks <- seq(0, 20, by = 5)
  tx <- x0 + bar_w * (ticks - VMIN) / (VMAX - VMIN)
  segments(ndc_x(tx), ndc_y(y0), ndc_x(tx), ndc_y(y0 - 0.006), col = TEXT_COLOR, lwd = 1)
  halo_text(ndc_x(tx), ndc_y(y0 - 0.010), ticks, upp, adj = c(0.5, 1), cex = 14 / 12)
  halo_text(ndc_x(x0 + w / 2), ndc_y(y0 - 0.048), "grams of carbon per square meter per week",
            upp, adj = c(0.5, 1), cex = 14 / 12)
  
  xc <- x0 + w / 2
  halo_text(ndc_x(xc), ndc_y(0.255), date_label, upp, adj = c(0.5, 0), cex = 22 / 12, font = 2)
  halo_text(ndc_x(xc), ndc_y(0.222), "Carbon Emissions", upp, adj = c(0.5, 0), cex = 18 / 12)
  text(ndc_x(0.99), ndc_y(0.012), CREDIT, adj = c(1, 0), col = "#bdbdbd", cex = 11 / 12)
}

render_frame <- function(path, L, bg_raster, fire, date) {
  open_png(path, L)
  rasterImage(bg_raster, L$x_left, L$y_bottom, L$x_right, L$y_top, interpolate = FALSE)
  rasterImage(fire, L$x_left, L$y_bottom, L$x_right, L$y_top, interpolate = FALSE)
  draw_legend(format(date, "%b  %Y"), L$mpp)
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
  
  L  <- map_layout()
  px <- pixel_lonlat(L)
  bg <- build_background(L, px)
  bg_raster <- as.raster(bg)
  lk  <- build_lookup(px, lat, lon)
  rm(px)
  lut <- make_lut()
  
  if (RENDER_MODE == "test") {
    dir.create(TEST_DIR, recursive = TRUE, showWarnings = FALSE)
    for (d in as.list(TEST_DATES)) {
      i <- which.min(abs(as.numeric(index$date - d)))
      nc <- nc_open(index$file[i])
      field <- read_week(nc, index$k[i])
      nc_close(nc)
      out <- file.path(TEST_DIR, sprintf("fires_foucaut_R_%s.png", index$date[i]))
      render_frame(out, L, bg_raster, fire_raster(field, lut, lk), index$date[i])
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
      render_frame(out, L, bg_raster, fire_raster(read_week(nc, index$k[i]), lut, lk), index$date[i])
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