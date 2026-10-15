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
# Script: saharan_dust_platecarree.R
# Description: Generates the Saharan dust optical depth animation from NASA
#     MERRA-2 hourly DUEXTTAU (2020-06-01 to 2020-06-30), rendered in the
#     PlateCarree projection.
# Note: For the data download, see download_merra2_dust_2020_06.py in the
#     saharan_dust scripts folder.
#
# Run:
#     # 1) Set PREVIEW_ONLY <- TRUE, check the PNG
#     # 2) Set PREVIEW_ONLY <- FALSE for the full render
#     caffeinate -i Rscript saharan_dust_platecarree.R
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
  library(rnaturalearth)
  library(av)
})
sf_use_s2(FALSE)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
VARIABLE  <- "DUEXTTAU"
DATA_PATH <- "/Volumes/CLIMATEDATA/dust/merra2_2020_06/MERRA2_DUEXTTAU_hourly_2020-06.nc"

LOCAL_DIR <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "saharan_dust", "_render")  # local APFS
FINAL_DIR <- "/Volumes/CLIMATEDATA/dust/animations"                                       # exFAT
NE_DIR    <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "saharan_dust", "_naturalearth")
OUT_NAME  <- "saharan_dust_june2020_platecarree_R.mp4"

PREVIEW_ONLY <- FALSE                     # TRUE: save one PNG and stop
PREVIEW_TIME <- "2020-06-06T10:30"        # matches the reference frame

# Map extent (lon_min, lon_max, lat_min, lat_max) — limited by the data subset
MAP_EXTENT <- c(-128, 42, 0, 55)

# Colorbar levels from the original (log-spaced, extended both ends)
LEVELS     <- c(0.10, 0.15, 0.21, 0.31, 0.45, 0.66, 0.97, 1.40, 2.10, 3.00)
CBAR_LABEL <- "Dust optical depth (550 nm)"

BG_COLOR     <- "#fdf5e8"                 # cream ocean/land background
CMAP_ANCHORS <- c("#fbe9cf", "#f6cf9c", "#eeac72", "#e2864f",
                  "#cd5f34", "#ad3c21", "#842312", "#5a0e06")
OVER_COLOR   <- "#3d0703"

SHOW_US_CITIES <- TRUE
CITY_MIN_POP   <- 50000

TITLE    <- "Saharan Dust, June 2020"
SUBTITLE <- "NASA MERRA-2 hourly dust optical depth (550 nm)"

LON_TICKS      <- seq(-120, 40, by = 20)  # 120°W ... 40°E
LAT_TICKS      <- seq(0, 50, by = 10)     # 0° ... 50°N
SHOW_GRIDLINES <- TRUE                    # FALSE: labels only, no lines

FPS     <- 24                             # 720 frames -> 30 s, same as original
FIGSIZE <- c(14, 5.2)                     # inches
DPI     <- 150                            # -> 2100 x 780 px
WIDTH   <- round(FIGSIZE[1] * DPI)
HEIGHT  <- round(FIGSIZE[2] * DPI)

# Figure-fraction positions (x0, x1, y0, y1), same layout as the Python version
AX_POS  <- c(0.045, 0.045 + 0.83, 0.08, 0.08 + 0.72)
CAX_POS <- c(0.915, 0.915 + 0.012, 0.12, 0.12 + 0.64)

PT2LWD <- 96 / 72                         # matplotlib points -> R lwd units

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
build_palette <- function() {
  n_bins <- length(LEVELS) - 1
  bins <- colorRampPalette(CMAP_ANCHORS)(n_bins)
  list(bins = bins, lookup = c(BG_COLOR, bins, OVER_COLOR))   # under, bins, over
}

# BoundaryNorm equivalent: < first level -> under, >= last level -> over, NA -> bad
values_to_raster <- function(z, lookup) {
  idx <- findInterval(z, LEVELS)                 # 0 = under, 1..9 bins, 10 = over
  cols <- lookup[idx + 1]
  cols[is.na(cols)] <- BG_COLOR
  m <- matrix(cols, nrow = nrow(z), ncol = ncol(z))   # [lon, lat], lat ascending
  as.raster(t(m)[ncol(z):1, , drop = FALSE])          # rows north -> south
}

parse_nc_time <- function(vals, units) {
  m <- regmatches(units, regexec("^\\s*(\\w+) since (.+)$", units))[[1]]
  mult <- c(seconds = 1, minutes = 60, hours = 3600, days = 86400)[[tolower(m[2])]]
  ref <- sub("(\\.\\d+)?\\s*(UTC|Z|[+-]\\d{2}:?\\d{2})?\\s*$", "", m[3])
  origin <- as.POSIXct(ref, tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S",
                                      "%Y-%m-%d %H:%M", "%Y-%m-%d"))
  origin + vals * mult
}

stamp <- function(t) format(t, "%Y-%m-%d %H:%Mz", tz = "UTC")

lon_label <- function(x) ifelse(x == 0, "0\u00b0",
                                paste0(abs(x), "\u00b0", ifelse(x < 0, "W", "E")))
lat_label <- function(x) ifelse(x == 0, "0\u00b0",
                                paste0(abs(x), "\u00b0", ifelse(x < 0, "S", "N")))

ne_layer <- function(scale, type, category) {
  dir.create(NE_DIR, recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    ne_load(scale = scale, type = type, category = category,
            destdir = NE_DIR, returnclass = "sf"),
    error = function(e)
      ne_download(scale = scale, type = type, category = category,
                  destdir = NE_DIR, load = TRUE, returnclass = "sf")
  )
}

crop_to_map <- function(x, pad = 2) {
  bb <- st_bbox(c(xmin = MAP_EXTENT[1] - pad, xmax = MAP_EXTENT[2] + pad,
                  ymin = MAP_EXTENT[3] - pad, ymax = MAP_EXTENT[4] + pad),
                crs = st_crs(4326))
  suppressWarnings(st_crop(st_make_valid(x), bb))
}

load_map_layers <- function() {
  layers <- list(
    coast   = crop_to_map(ne_layer(50, "coastline", "physical")),
    borders = crop_to_map(ne_layer(50, "admin_0_boundary_lines_land", "cultural")),
    states  = crop_to_map(ne_layer(50, "admin_1_states_provinces_lakes", "cultural")),
    cities  = NULL
  )
  if (SHOW_US_CITIES) {
    layers$cities <- tryCatch({
      pp <- ne_layer(10, "populated_places", "cultural")
      names(pp)[names(pp) != attr(pp, "sf_column")] <-
        toupper(names(pp)[names(pp) != attr(pp, "sf_column")])
      pp <- pp[pp$ADM0_A3 == "USA" & !is.na(pp$POP_MAX) & pp$POP_MAX >= CITY_MIN_POP, ]
      st_coordinates(pp)
    }, error = function(e) {
      message("City layer skipped (", conditionMessage(e), ")")
      NULL
    })
  }
  layers
}

# Overlay a full-figure [0,1] x [0,1] coordinate system for fig.text equivalents
fig_text <- function(x, y, label, adj, cex, font = 1, col = "black") {
  par(plt = c(0, 1, 0, 1), new = TRUE)
  plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  text(x, y, label, adj = adj, cex = cex, font = font, col = col, xpd = NA)
}

draw_frame <- function(img, img_ext, layers, pal, stamp_text) {
  # --- Map axes -------------------------------------------------------------
  par(plt = AX_POS, xaxs = "i", yaxs = "i", mgp = c(3, 0.3, 0))
  plot.new()
  plot.window(xlim = MAP_EXTENT[1:2], ylim = MAP_EXTENT[3:4])
  rect(MAP_EXTENT[1], MAP_EXTENT[3], MAP_EXTENT[2], MAP_EXTENT[4],
       col = BG_COLOR, border = NA)
  rasterImage(img, img_ext[1], img_ext[3], img_ext[2], img_ext[4], interpolate = TRUE)
  
  plot(st_geometry(layers$coast),   add = TRUE, col = "#555555", lwd = 0.4 * PT2LWD)
  plot(st_geometry(layers$borders), add = TRUE, col = "#777777", lwd = 0.3 * PT2LWD)
  plot(st_geometry(layers$states),  add = TRUE, col = NA, border = "#777777",
       lwd = 0.3 * PT2LWD)
  
  if (SHOW_GRIDLINES) {
    abline(v = LON_TICKS, h = LAT_TICKS, col = adjustcolor("#999999", 0.6),
           lty = "dashed", lwd = 0.3 * PT2LWD)
  }
  if (!is.null(layers$cities) && nrow(layers$cities) > 0) {
    points(layers$cities[, 1], layers$cities[, 2], pch = 16, cex = 0.2, col = "black")
  }
  
  axis(1, at = LON_TICKS, labels = lon_label(LON_TICKS), tick = FALSE,
       cex.axis = 11 / 12, col.axis = "#333333")
  axis(2, at = LAT_TICKS, labels = lat_label(LAT_TICKS), tick = FALSE, las = 1,
       cex.axis = 11 / 12, col.axis = "#333333")
  
  # --- Colorbar (uniform bins, triangular extensions) -----------------------
  nb <- length(pal$bins)
  ext <- 0.05 * nb                                   # matplotlib extendfrac = 0.05
  par(plt = CAX_POS, new = TRUE, mgp = c(3, 0.4, 0))
  plot.new()
  plot.window(c(0, 1), c(0, nb), xaxs = "i", yaxs = "i")
  rect(0, 0:(nb - 1), 1, 1:nb, col = pal$bins, border = NA)
  polygon(c(0, 0.5, 1), c(0, -ext, 0), col = BG_COLOR, border = NA, xpd = NA)
  polygon(c(0, 0.5, 1), c(nb, nb + ext, nb), col = OVER_COLOR, border = NA, xpd = NA)
  polygon(c(0, 0, 0.5, 1, 1, 0.5), c(0, nb, nb + ext, nb, 0, -ext),
          col = NA, border = "black", lwd = 0.8 * PT2LWD, xpd = NA)
  axis(4, at = 0:nb, labels = sprintf("%.2f", LEVELS), las = 1, lwd = 0,
       lwd.ticks = 0.8 * PT2LWD, tcl = -0.25, cex.axis = 11 / 12)
  mtext(CBAR_LABEL, side = 2, line = 0.5, cex = 12 / 12)
  
  # --- Figure text ----------------------------------------------------------
  fig_text(0.46, 0.925, TITLE, adj = c(0.5, 0), cex = 15 / 12, font = 2)
  fig_text(0.46, 0.875, SUBTITLE, adj = c(0.5, 0), cex = 10 / 12, col = "#444444")
  fig_text(0.985, 0.025, stamp_text, adj = c(1, 0), cex = 9 / 12)
}

open_png <- function(path) {
  png(path, width = WIDTH, height = HEIGHT, res = DPI, pointsize = 12, bg = "white")
}

verify_mp4 <- function(path, expected_frames) {
  if (!file.exists(path) || file.size(path) == 0) stop(basename(path), " is missing or empty")
  info <- av_media_info(path)
  n_frames <- info$video$frames
  if (abs(n_frames - expected_frames) > 1) {
    stop(sprintf("%s: %d frames, expected %d", basename(path), n_frames, expected_frames))
  }
  cat(sprintf("Verified %s: ~%d frames, %.2f s (expected %.2f s)\n",
              basename(path), n_frames, info$duration, expected_frames / FPS))
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  if (!file.exists(DATA_PATH)) {
    stop(DATA_PATH, " not found. Is /Volumes/CLIMATEDATA mounted?")
  }
  dir.create(LOCAL_DIR, recursive = TRUE, showWarnings = FALSE)
  
  nc <- nc_open(DATA_PATH)
  on.exit(nc_close(nc), add = TRUE)
  
  var <- nc$var[[VARIABLE]]
  dim_names <- vapply(var$dim, function(d) d$name, character(1))
  lon_all <- ncvar_get(nc, "lon")
  lat_all <- ncvar_get(nc, "lat")
  times <- parse_nc_time(ncvar_get(nc, "time"), ncatt_get(nc, "time", "units")$value)
  
  # Spatial subset (same as the Python .sel slices)
  i_lon <- which(lon_all >= MAP_EXTENT[1] - 1 & lon_all <= MAP_EXTENT[2] + 1)
  i_lat <- which(lat_all >= MAP_EXTENT[3] & lat_all <= MAP_EXTENT[4] + 1)
  lon <- lon_all[i_lon]
  lat <- lat_all[i_lat]
  n <- length(times)
  cat(sprintf("Loaded %s: %d hourly frames, %d x %d\n", VARIABLE, n, length(lat), length(lon)))
  
  lat_desc <- lat[1] > lat[length(lat)]
  if (lat_desc) lat <- rev(lat)
  
  read_frame <- function(k) {
    start <- c(lon = min(i_lon), lat = min(i_lat), time = k)[dim_names]
    count <- c(lon = length(i_lon), lat = length(i_lat), time = 1)[dim_names]
    z <- ncvar_get(nc, var, start = start, count = count, collapse_degen = TRUE)
    if (which(dim_names == "lat") < which(dim_names == "lon")) z <- t(z)  # -> [lon, lat]
    if (lat_desc) z <- z[, ncol(z):1, drop = FALSE]
    z
  }
  
  dlon <- lon[2] - lon[1]
  dlat <- lat[2] - lat[1]
  img_ext <- c(lon[1] - dlon / 2, lon[length(lon)] + dlon / 2,
               lat[1] - dlat / 2, lat[length(lat)] + dlat / 2)
  
  pal <- build_palette()
  layers <- load_map_layers()
  
  if (PREVIEW_ONLY) {
    target <- as.POSIXct(PREVIEW_TIME, format = "%Y-%m-%dT%H:%M", tz = "UTC")
    k <- which.min(abs(as.numeric(difftime(times, target, units = "secs"))))
    out <- file.path(LOCAL_DIR, "preview_frame_R.png")
    open_png(out)
    draw_frame(values_to_raster(read_frame(k), pal$lookup), img_ext, layers, pal,
               stamp(times[k]))
    dev.off()
    cat("Preview saved:", out, "\n")
    return(invisible())
  }
  
  # Render frames to local APFS
  frame_dir <- file.path(LOCAL_DIR, "frames_R")
  unlink(frame_dir, recursive = TRUE)
  dir.create(frame_dir, recursive = TRUE)
  frame_files <- file.path(frame_dir, sprintf("frame_%04d.png", seq_len(n)))
  
  for (k in seq_len(n)) {
    open_png(frame_files[k])
    draw_frame(values_to_raster(read_frame(k), pal$lookup), img_ext, layers, pal,
               stamp(times[k]))
    dev.off()
    if ((k - 1) %% 24 == 0) cat(sprintf("  frame %d/%d  %s\n", k, n, stamp(times[k])))
  }
  
  # Encode, verify, then move to the exFAT drive
  local_mp4 <- file.path(LOCAL_DIR, OUT_NAME)
  if (file.exists(local_mp4)) file.remove(local_mp4)
  av_encode_video(frame_files, output = local_mp4, framerate = FPS,
                  codec = "libx264", vfilter = "format=yuv420p", verbose = FALSE)
  verify_mp4(local_mp4, n)
  
  dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
  final_mp4 <- file.path(FINAL_DIR, OUT_NAME)
  if (file.exists(final_mp4)) file.remove(final_mp4)
  if (!file.copy(local_mp4, final_mp4) || file.size(final_mp4) != file.size(local_mp4)) {
    stop("Copy to ", FINAL_DIR, " failed; local file kept at ", local_mp4)
  }
  file.remove(local_mp4)
  unlink(frame_dir, recursive = TRUE)
  cat("Saved", final_mp4, "\n")
}

main()