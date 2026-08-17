# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: north_stereo.R
# Description: Generates the Northern Hemisphere seasonal-cycle stereographic animation
# from NASA Blue Marble Next Generation monthly imagery (2004) with an ERA5 sea_ice_cover
# overlay (24 biweekly frames, 1st and 15th of each month), rendered in the North Polar
# Stereographic projection.
# Note: For the Southern Hemisphere, see south_stereo.R in the stereographic scripts folder.

suppressPackageStartupMessages({
  library(terra)
  library(av)
})

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

BLUE_MARBLE_DIR <- "/Volumes/CLIMATEDATA/blue_marble/2004"
ICE_FILE <- "/Volumes/CLIMATEDATA/era5_sea_ice/2004/era5_sea_ice_cover_2004_biweekly.nc"

# R outputs are isolated in r_output/ to avoid overwriting Python-generated files
FRAMES_DIR <- "/Volumes/CLIMATEDATA/stereographic/north/r_output/frames"
VIDEO_OUT <- "/Volumes/CLIMATEDATA/stereographic/north/r_output/north_seasonal_cycle_r.mp4"

TARGET_DATES <- as.Date(unlist(lapply(1:12, function(m) {
  sprintf("2004-%02d-%02d", m, c(1, 15))
})))

ICE_THRESHOLD <- 0.15   # minimum concentration to render as visible ice edge
NORTH_MIN_LAT <- 0      # southern boundary of the hemisphere view (equator)
NORTH_CENTRAL_LON <- 300  # meridian pointing to the top of the map. PROJ's stereographic
# +lon_0 does not necessarily rotate the same direction as
# Cartopy's central_longitude -- this value was carried over
# from the tuned Python script as a starting guess and should
# be retuned by eye rather than assumed correct.
FPS <- 6

CRS_NORTH <- sprintf("+proj=stere +lat_0=90 +lon_0=%d +datum=WGS84 +units=m +no_defs",
                     NORTH_CENTRAL_LON)

# Template grid the base image and ice layer are both projected onto, so the
# two layers land on identical pixel grids. ~13,000 km half-width comfortably
# covers the full hemisphere disc (equator sits ~12,740 km from the pole in
# true polar stereographic); 15 km resolution balances detail vs. render time.
TEMPLATE_NORTH <- terra::rast(
  xmin = -13000000, xmax = 13000000,
  ymin = -13000000, ymax = 13000000,
  resolution = 15000,
  crs = CRS_NORTH
)

# ------------------------------------------------------------------------------
# Blue Marble loading and biweekly blending
# ------------------------------------------------------------------------------

month_file <- function(month) {
  file.path(BLUE_MARBLE_DIR, sprintf("world.2004%02d.3x5400x2700.jpg", month))
}

load_month_image <- function(month) {
  img <- terra::rast(month_file(month))
  terra::ext(img) <- terra::ext(-180, 180, -90, 90)
  terra::crs(img) <- "EPSG:4326"
  # Flip applied after ext/crs assignment -- doing it before appeared to have
  # no effect, possibly because reassigning ext() on a lazy raster discarded
  # the pending flip.
  img <- terra::flip(img, direction = "vertical")
  img
}

blend_base_image <- function(target_date) {
  # Blue Marble is monthly, not biweekly. Blend the two nearest monthly
  # composites (each anchored to the 15th of its month) to build a base
  # image for the target biweekly date. Dates outside the first/last
  # anchor fall back to the nearest single month with no blending.
  centers <- as.Date(sprintf("2004-%02d-15", 1:12))
  names(centers) <- 1:12
  
  if (target_date <= centers[["1"]]) return(load_month_image(1))
  if (target_date >= centers[["12"]]) return(load_month_image(12))
  
  prev_month <- max(as.integer(names(centers)[centers <= target_date]))
  next_month <- prev_month + 1
  prev_center <- centers[[as.character(prev_month)]]
  next_center <- centers[[as.character(next_month)]]
  
  if (target_date == prev_center) return(load_month_image(prev_month))
  
  span <- as.numeric(next_center - prev_center)
  weight <- as.numeric(target_date - prev_center) / span
  
  prev_img <- load_month_image(prev_month)
  next_img <- load_month_image(next_month)
  (1 - weight) * prev_img + weight * next_img
}

# ------------------------------------------------------------------------------
# ERA5 sea ice loading
# ------------------------------------------------------------------------------

load_ice_dataset <- function() {
  r <- terra::rast(ICE_FILE)
  terra::rotate(r)  # normalize ERA5's 0-360 longitude to -180/180
}

get_ice_frame <- function(ice_stack, target_date) {
  times <- terra::time(ice_stack)
  idx <- which.min(abs(as.numeric(as.Date(times) - target_date)))
  ice_stack[[idx]]
}

# ------------------------------------------------------------------------------
# Rendering
# ------------------------------------------------------------------------------

hemisphere_crop <- function(r) {
  terra::crop(r, terra::ext(-180, 180, NORTH_MIN_LAT, 90))
}

render_frame <- function(target_date, base_img, ice_layer, output_path) {
  base_cropped <- hemisphere_crop(base_img)
  ice_cropped <- hemisphere_crop(ice_layer)
  
  base_proj <- terra::project(base_cropped, TEMPLATE_NORTH, method = "bilinear")
  ice_proj <- terra::project(ice_cropped, TEMPLATE_NORTH, method = "near")
  
  # NA below threshold (fully transparent); alpha ramps with concentration above it
  ice_masked <- terra::classify(ice_proj, rbind(c(-Inf, ICE_THRESHOLD, NA)))
  alpha_ramp <- grDevices::colorRampPalette(
    c(grDevices::rgb(1, 1, 1, 0.3), grDevices::rgb(1, 1, 1, 1)), alpha = TRUE
  )(100)
  
  grDevices::png(output_path, width = 1200, height = 1200, res = 150, bg = "white")
  par(bg = "white", mar = c(4, 2, 6, 2))
  
  terra::plotRGB(base_proj, r = 1, g = 2, b = 3, scale = 255,
                 axes = FALSE, mar = c(4, 2, 6, 2), colNA = "white")
  terra::plot(ice_masked, col = alpha_ramp, range = c(ICE_THRESHOLD, 1),
              add = TRUE, legend = FALSE, axes = FALSE)
  
  title(main = "Northern Hemisphere Seasonal Cycle", cex.main = 1.6, line = 3.2)
  mtext("NASA Blue Marble (2004) base imagery with sea ice from ERA5 reanalysis",
        side = 3, line = 1, cex = 0.7)
  mtext(format(target_date, "%b %d, %Y"), side = 1, line = 1, adj = 1, cex = 0.9)
  
  grDevices::dev.off()
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

main <- function() {
  dir.create(FRAMES_DIR, recursive = TRUE, showWarnings = FALSE)
  ice_ds <- load_ice_dataset()
  
  frame_paths <- character(length(TARGET_DATES))
  for (i in seq_along(TARGET_DATES)) {
    target_date <- TARGET_DATES[i]
    cat(sprintf("Rendering frame: %s\n", format(target_date, "%Y-%m-%d")))
    
    base_img <- blend_base_image(target_date)
    ice_layer <- get_ice_frame(ice_ds, target_date)
    
    out_path <- file.path(FRAMES_DIR, sprintf("north_%s.png", format(target_date, "%Y-%m-%d")))
    render_frame(target_date, base_img, ice_layer, out_path)
    frame_paths[i] <- out_path
  }
  
  cat(sprintf("Assembling %d frames into %s...\n", length(frame_paths), VIDEO_OUT))
  dir.create(dirname(VIDEO_OUT), recursive = TRUE, showWarnings = FALSE)
  av::av_encode_video(frame_paths, output = VIDEO_OUT, framerate = FPS)
  
  cat("Done.\n")
}

main()