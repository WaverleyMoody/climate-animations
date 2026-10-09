# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations
# Library by Professor John Michael Wallace.
#
# Script: animate_fires_platecarree.R
#
# Description: Generates the fire carbon emissions time-series animation from
#              GFED5.1 weekly totals (March 2003 - January 2022), rendered in the
#              Plate Carree projection in the style of the NASA SVS "Fires - a
#              global perspective" animation.
#
# Note: For the Robinson, Foucaut, and Nicolosi projections, see the other scripts
#       in the fires scripts folder.
#
# Run:
#   1. RENDER_MODE <- "test"  -> writes a few PNG frames for review (fast)
#   2. RENDER_MODE <- "full"  -> caffeinate -i Rscript animate_fires_platecarree.R
# Frames are written to FRAME_DIR and skipped if already present, so an
# interrupted full render resumes where it stopped.

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
FRAME_DIR   <- "/Volumes/CLIMATEDATA/_frames_tmp_fires_platecarree_R"
BG_CACHE    <- file.path(WORK, "fires_background_platecarree_R.png")
FINAL_VIDEO <- "/Volumes/CLIMATEDATA/fires_2003_2022_platecarree_R.mp4"

WIDTH <- 1920; HEIGHT <- 1080; RES <- 100   # px; RES makes font sizes match the Python version
FPS   <- 12                                 # 987 weeks -> ~82 s (NASA original runs 1:24)

# Basemap colors (NASA SVS look)
OCEAN_COLOR  <- "#4d5868"
LAND_DARK    <- 0.10; LAND_LIGHT <- 0.36    # gray range for shaded relief (0 = black, 1 = white)
LAND_FLAT    <- "#2e2e2e"                   # used only if the relief raster is unavailable
BORDER_COLOR <- adjustcolor("#d9d9d9", alpha.f = 0.6)
RELIEF_URLS  <- c(
  "https://naciscdn.org/naturalearth/50m/raster/GRAY_50M_SR.zip",
  "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/GRAY_50M_SR.zip"
)

# Data colors: green -> orange -> white, 0-20 g C m-2 week-1, saturating above 20
VMIN <- 0; VMAX <- 20
CMAP_NODES <- c("#5e8c3a", "#a3b257", "#e5b260", "#f0d5a4", "#f3f1ec")   # at 0, .25, .5, .75, 1
ALPHA_MIN  <- 0.05; ALPHA_FULL <- 1.0      # transparent below ALPHA_MIN, opaque above ALPHA_FULL
TEXT_COLOR <- "#e6e6e6"
CREDIT     <- "Data: GFED5.1 (van der Werf et al.)  |  SDSU Climate Informatics Lab"

LUT_STEP <- 0.02                            # g C m-2 week-1 per color-table bin

invisible(Sys.setlocale("LC_TIME", "C"))    # English month abbreviations
PNG_TYPE <- if (capabilities("cairo")) "cairo" else NULL


# ---------------------------------------------------------------- helpers
move_file <- function(from, to) {
  # Cross-filesystem move (local APFS -> exFAT): copy, check size, then delete
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok || file.size(to) != file.size(from)) stop("Copy to ", to, " failed; left at ", from)
  unlink(from)
}

open_png <- function(path) {
  args <- list(filename = path, width = WIDTH, height = HEIGHT, res = RES, bg = OCEAN_COLOR)
  if (!is.null(PNG_TYPE)) args$type <- PNG_TYPE
  do.call(png, args)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(-180, 180), ylim = c(-90, 90))   # no asp: stretched to 16:9 like NASA
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
  # Natural Earth 1:50m gray shaded relief (shared with the Python scripts); NULL if unavailable
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
  r <- suppressWarnings(rast(tif))[[1]]               # the .tif has no embedded georeference
  ext(r) <- ext(-180, 180, -90, 90)
  crs(r) <- "EPSG:4326"
  # Average down to exactly one value per output pixel
  template <- rast(ncols = WIDTH, nrows = HEIGHT, xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                   crs = "EPSG:4326")
  resample(r, template, method = "average")
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

build_background <- function() {
  # Render the static basemap once at output size; returns a (HEIGHT, WIDTH, 3) array in 0-1
  if (file.exists(BG_CACHE)) return(readPNG(BG_CACHE)[, , 1:3])
  
  message("Building background ...")
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  ocean <- st_geometry(ne_layer("ocean", "physical"))
  lakes <- st_geometry(ne_layer("lakes", "physical"))
  borders <- st_geometry(ne_layer("admin_0_boundary_lines_land", "cultural"))
  relief <- load_relief()
  
  open_png(BG_CACHE)
  if (!is.null(relief)) {
    shade <- LAND_DARK + (LAND_LIGHT - LAND_DARK) * (as.matrix(relief, wide = TRUE) / 255)
    shade[is.na(shade)] <- LAND_DARK
    rasterImage(as.raster(matrix(gray(shade), nrow = nrow(shade))), -180, -90, 180, 90,
                interpolate = FALSE)
  } else {
    land <- st_geometry(ne_layer("land", "physical"))
    plot(land, col = LAND_FLAT, border = NA, add = TRUE)
  }
  plot(ocean, col = OCEAN_COLOR, border = NA, add = TRUE)
  plot(lakes, col = OCEAN_COLOR, border = NA, add = TRUE)
  plot(borders, col = BORDER_COLOR, lwd = 0.9, add = TRUE)
  dev.off()
  
  bg <- readPNG(BG_CACHE)[, , 1:3]
  if (any(dim(bg)[1:2] != c(HEIGHT, WIDTH))) stop("Background has wrong size: ", paste(dim(bg), collapse = "x"))
  bg
}


# ---------------------------------------------------------------- colors and legend
make_lut <- function() {
  # One color (with alpha) per LUT_STEP bin from 0 to VMAX; the last bin is the >= VMAX color
  v <- seq(0, VMAX, by = LUT_STEP)
  ramp <- colorRamp(CMAP_NODES)
  rgb_v <- ramp(pmin(pmax((v - VMIN) / (VMAX - VMIN), 0), 1))
  alpha <- pmin(pmax((v - ALPHA_MIN) / (ALPHA_FULL - ALPHA_MIN), 0), 1)
  rgb(rgb_v[, 1], rgb_v[, 2], rgb_v[, 3], alpha * 255, maxColorValue = 255)
}

fire_raster <- function(field, lut, lat_ascending) {
  # field: [lon, lat] matrix from ncdf4 -> north-up raster of hex colors
  idx <- pmin(floor(field / LUT_STEP) + 1, length(lut))
  idx[!is.finite(idx) | idx < 1] <- 1
  m <- matrix(lut[idx], nrow = nrow(field))         # [lon, lat]
  m <- t(m)                                          # [lat, lon]
  if (lat_ascending) m <- m[nrow(m):1, ]             # north at the top
  as.raster(m)
}

ndc_x <- function(f) grconvertX(f, "ndc", "user")
ndc_y <- function(f) grconvertY(f, "ndc", "user")

draw_legend <- function(date_label) {
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
  text(ndc_x(tx), ndc_y(y0 - 0.010), ticks, adj = c(0.5, 1), col = TEXT_COLOR, cex = 14 / 12)
  text(ndc_x(x0 + w / 2), ndc_y(y0 - 0.048), "grams of carbon per square meter per week",
       adj = c(0.5, 1), col = TEXT_COLOR, cex = 14 / 12)
  
  xc <- x0 + w / 2
  text(ndc_x(xc), ndc_y(0.255), date_label, adj = c(0.5, 0), col = TEXT_COLOR, cex = 22 / 12, font = 2)
  text(ndc_x(xc), ndc_y(0.222), "Carbon Emissions", adj = c(0.5, 0), col = TEXT_COLOR, cex = 18 / 12)
  text(ndc_x(0.99), ndc_y(0.012), CREDIT, adj = c(1, 0), col = "#bdbdbd", cex = 11 / 12)
}

render_frame <- function(path, bg_raster, fire, date) {
  open_png(path)
  rasterImage(bg_raster, -180, -90, 180, 90, interpolate = FALSE)
  rasterImage(fire, -180, -90, 180, 90, interpolate = FALSE)
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
  lat <- ncvar_get(nc0, "lat")
  nc_close(nc0)
  lat_ascending <- lat[1] < lat[length(lat)]
  
  bg <- build_background()
  bg_raster <- as.raster(bg)
  lut <- make_lut()
  
  if (RENDER_MODE == "test") {
    dir.create(TEST_DIR, recursive = TRUE, showWarnings = FALSE)
    for (d in as.list(TEST_DATES)) {
      i <- which.min(abs(as.numeric(index$date - d)))
      nc <- nc_open(index$file[i])
      field <- read_week(nc, index$k[i])
      nc_close(nc)
      out <- file.path(TEST_DIR, sprintf("fires_platecarree_R_%s.png", index$date[i]))
      render_frame(out, bg_raster, fire_raster(field, lut, lat_ascending), index$date[i])
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
      render_frame(out, bg_raster, fire_raster(read_week(nc, index$k[i]), lut, lat_ascending),
                   index$date[i])
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