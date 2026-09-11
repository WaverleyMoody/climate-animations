# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation Animations Library
# by Professor John Michael Wallace.
#
# Script: animate_sst_orthographic.R

# Description: Generates the sea surface temperature (SST) daily time-series animation from
#   JPL MUR25 (2016-2020), rendered in the orthographic projection as a spinning globe that
#   gradually tilts to reveal the South Pole over the course of the animation.

# Note: Translated from animate_sst_orthographic.py. 

library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(grid)
library(tidyterra)
library(av)
library(glue)
library(progress)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
IN_PATH  <- "/Volumes/CLIMATEDATA/sst_mur25/timeseries/sst_daily_timeseries_2016_2020_025deg_celsius.nc"
OUT_DIR  <- "/Users/waverleymoody/CLIMATE ANIMATIONS/sea_surface_temp/_local_output"
OUT_PATH <- file.path(OUT_DIR, "sst_daily_2016_2020_orthographic.mp4")

# Frame-by-frame PNGs are written here, then assembled into OUT_PATH -- R has
# no single call that both renders and encodes in one pass the way
# matplotlib's FuncAnimation + FFMpegWriter does, so this is a two-stage
# process bridged by av::av_encode_video() at the end.
FRAMES_DIR <- file.path(OUT_DIR, "_frames")

VAR_NAME <- "analysed_sst"

VMIN <- -2; VMAX <- 32          # degrees C
CMAP_COLORS <- rev(c(           # approximates matplotlib's "RdYlBu_r"
  "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090",
  "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"
))

START_LAT <- 20                 # initial viewing latitude (tilt of the globe)
END_LAT <- -55                  # final viewing latitude, tipped toward South Pole
TILT_START_FRAC <- 0.6          # tilt spans the whole animation for a gradual reveal

ROTATIONS <- 2                  # number of full 360-degree spins over the whole animation

FIG_WIDTH_PX <- 1200; FIG_HEIGHT_PX <- 1200
DPI <- 150

TARGET_RUNTIME_SEC <- 104       # 1:44

TITLE_TEXT <- "Sea Surface Temperature"

# Set to a frame index (0-indexed, e.g. 0 for the very first frame) to render
# ONLY that one frame and stop -- skips the full loop and the video encoding
# step entirely. Set back to NULL to run the full animation.
TEST_FRAME_INDEX <- NULL
# ---------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FRAMES_DIR, showWarnings = FALSE, recursive = TRUE)


.verify_readable <- function(path) {
  result <- tryCatch({
    r <- terra::rast(path, subds = VAR_NAME)
    terra::values(r[[1]])  # force a real read of the first layer, not just metadata
    TRUE
  }, error = function(e) {
    message(glue::glue("Verification failed: {conditionMessage(e)}"))
    FALSE
  })
  result
}


.format_date <- function(posix_time) {
  # Equivalent to the Python version's manual "%b %-d, %Y" workaround --
  # format() with "%e" (space-padded day) plus trimws() avoids the same
  # platform-specific leading-zero/no-leading-zero inconsistency.
  day <- trimws(format(posix_time, "%e"))
  glue::glue("{format(posix_time, '%b')} {day}, {format(posix_time, '%Y')}")
}


.visible_hemisphere_mask <- function(lon_i, lat_i) {
  # sf::st_transform() has no concept of "the far side of the globe" -- it
  # will happily reproject a polygon that straddles the horizon (90 degrees
  # from the view center) into garbled, self-intersecting coordinates, which
  # ggplot then silently fails to draw. This produced the missing USA /
  # Eastern Europe polygons.
  #
  # Fix: build the visible hemisphere as an actual polygon in lon/lat space
  # by drawing a circle of Earth's radius around the origin IN the ortho
  # projection (which is exactly the horizon), then transforming that circle
  # back to lon/lat. Intersecting land with this polygon BEFORE reprojecting
  # forward clips every country to what's actually visible, so no polygon
  # ever crosses the horizon in the first place.
  ortho_crs <- glue::glue(
    "+proj=ortho +lat_0={lat_i} +lon_0={lon_i} +datum=WGS84 +units=m +no_defs"
  )
  earth_radius_m <- 6371000
  circle_ortho <- sf::st_buffer(
    sf::st_sfc(sf::st_point(c(0, 0)), crs = ortho_crs),
    dist = earth_radius_m - 1000  # slightly inside the true horizon to avoid edge-case numerical issues
  )
  sf::st_transform(circle_ortho, 4326)
}


render_one_frame <- function(frame_idx0, sst, times, land, central_lons, central_lats) {
  # frame_idx0 is 0-indexed (matching the Python frame_idx); R vectors/lists
  # below are 1-indexed, hence the +1 throughout.
  i <- frame_idx0 + 1
  frame_path <- file.path(FRAMES_DIR, sprintf("frame_%05d.png", frame_idx0))
  
  lon_i <- central_lons[i]
  lat_i <- central_lats[i]
  ortho_crs <- glue::glue(
    "+proj=ortho +lat_0={lat_i} +lon_0={lon_i} +datum=WGS84 +units=m +no_defs"
  )
  
  # Reproject this frame's SST layer and the land polygons into the
  # current orthographic view. terra::project() naturally leaves the
  # far hemisphere as NA, giving the same "hidden back side" effect
  # cartopy's Orthographic produces automatically.
  sst_layer <- sst[[i]]
  sst_ortho <- terra::project(sst_layer, ortho_crs)
  
  land_ortho <- tryCatch({
    hemisphere_mask <- .visible_hemisphere_mask(lon_i, lat_i)
    land_visible <- suppressWarnings(
      sf::st_intersection(sf::st_make_valid(land), hemisphere_mask)
    )
    sf::st_transform(land_visible, ortho_crs)
  }, error = function(e) NULL)  # if clipping itself fails for some reason, skip land for this frame
  
  p <- ggplot() +
    tidyterra::geom_spatraster(data = sst_ortho) +
    scale_fill_gradientn(
      colors = CMAP_COLORS, limits = c(VMIN, VMAX),
      oob = scales::squish, na.value = "black",
      name = "Sea Surface Temperature (\u00b0C)",
      guide = guide_colorbar(
        barwidth = unit(0.75 * FIG_WIDTH_PX / DPI, "in"),  # ~75% of figure width, much longer than the default
        barheight = unit(0.12, "in"),
        title.position = "top"
      )
    )
  
  if (!is.null(land_ortho)) {
    p <- p + geom_sf(data = land_ortho, fill = "#D3D3D3", color = "black", linewidth = 0.2)
  }
  
  p <- p +
    coord_sf(crs = ortho_crs, datum = NA) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.position = "bottom",
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", face = "bold", hjust = 0.03, size = 13),
      plot.subtitle = element_text(color = "white", face = "bold", hjust = 0.97, size = 13),
    ) +
    labs(title = TITLE_TEXT, subtitle = .format_date(times[i]))
  
  ggsave(
    frame_path, plot = p,
    width = FIG_WIDTH_PX / DPI, height = FIG_HEIGHT_PX / DPI,
    dpi = DPI, bg = "black"
  )
  
  frame_path
}


main <- function() {
  message("Opening dataset (lazy per-layer access via terra)...")
  sst <- terra::rast(IN_PATH, subds = VAR_NAME)
  
  n_frames <- terra::nlyr(sst)
  times <- terra::time(sst)  # POSIXct vector, one per layer
  
  fps <- n_frames / TARGET_RUNTIME_SEC
  message(glue::glue("Total frames: {n_frames}"))
  message(glue::glue("FPS set to {round(fps, 2)} to hit target runtime of {TARGET_RUNTIME_SEC}s"))
  
  # ------------------------------------------------------------------
  # Precompute per-frame longitude (continuous spin) and latitude (tilt),
  # matching the Python version's precomputed central_lons/central_lats.
  # ------------------------------------------------------------------
  lon_per_frame <- (ROTATIONS * 360.0) / n_frames
  central_lons <- sapply(0:(n_frames - 1), function(i) {
    ((-i * lon_per_frame) %% 360)
  })
  
  tilt_start_frame <- as.integer(n_frames * TILT_START_FRAC)
  central_lats <- sapply(0:(n_frames - 1), function(i) {
    if (i < tilt_start_frame) {
      START_LAT
    } else {
      progress_frac <- (i - tilt_start_frame) / max(1, (n_frames - 1 - tilt_start_frame))
      eased <- 0.5 - 0.5 * cos(pi * progress_frac)  # smoothstep easing, same as Python version
      START_LAT + eased * (END_LAT - START_LAT)
    }
  })
  
  # Land polygons, fetched once and reprojected per frame below (mirrors
  # cartopy's cfeature.LAND, which is likewise a static vector layer
  # transformed into whatever the current GeoAxes projection is).
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # --- Single-frame test mode -----------------------------------------
  if (!is.null(TEST_FRAME_INDEX)) {
    message(glue::glue("TEST MODE: rendering only frame {TEST_FRAME_INDEX}, then stopping."))
    if (TEST_FRAME_INDEX < 0 || TEST_FRAME_INDEX >= n_frames) {
      stop(glue::glue("TEST_FRAME_INDEX must be between 0 and {n_frames - 1}"))
    }
    result_path <- render_one_frame(
      frame_idx0 = TEST_FRAME_INDEX,
      sst = sst, times = times, land = land,
      central_lons = central_lons, central_lats = central_lats
    )
    message(glue::glue("Test frame saved to: {result_path}"))
    message("Open that file to check colors, projection, land overlay, and title/date text.")
    return(invisible(NULL))
  }
  # ---------------------------------------------------------------------
  
  pb <- progress::progress_bar$new(
    format = "Rendering frame :current/:total [:bar] :percent eta: :eta",
    total = n_frames, clear = FALSE, width = 80
  )
  
  for (i in seq_len(n_frames)) {
    frame_idx0 <- i - 1  # 0-indexed, matching the Python frame_idx
    render_one_frame(
      frame_idx0 = frame_idx0,
      sst = sst, times = times, land = land,
      central_lons = central_lons, central_lats = central_lats
    )
    pb$tick()
  }
  
  message(glue::glue("Encoding {n_frames} frames to {OUT_PATH} at {round(fps, 2)} fps..."))
  frame_files <- sprintf(file.path(FRAMES_DIR, "frame_%05d.png"), 0:(n_frames - 1))
  av::av_encode_video(frame_files, output = OUT_PATH, framerate = fps)
  
  message("Verifying output file...")
  if (file.exists(OUT_PATH) && file.info(OUT_PATH)$size > 0) {
    message(glue::glue("Done. Output saved: {OUT_PATH}"))
  } else {
    warning("Output file missing or empty -- check for errors above.")
  }
}

main()