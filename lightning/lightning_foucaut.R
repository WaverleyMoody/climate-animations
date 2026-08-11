# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: lightning_foucaut.R
# Description: Generates the lightning climatology animation from WWLLN/WGLC data (2010-2025), rendered in the Foucaut projection. 365 calendar-day frames (Feb 29 excluded), each the multi-year daily-climatological mean, annualized (x365) to match the reference product's strokes km-2 yr-1 units.
#
# Note: For the Plate Carrée, Robinson, and Nicolosi projections, see the other scripts in the lightning scripts folder.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Volumes/CLIMATEDATA/lightning/climatology/climatology_lightning_365frame.nc"
output_dir <- "/Users/waverleymoody/Downloads/lightning_animation/r_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

proj_name  <- "foucaut"
target_crs <- "+proj=fouc +datum=WGS84 +units=m +no_defs"

frames_dir <- file.path(output_dir, proj_name, "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

# ── Plot settings (log scale) ────────────────────────────────────────────
vmin <- 0.003
vmax <- 30
tick_values <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)

# log_pos(value) maps a data value to its fractional position (0-1) along
# the log10 range from vmin to vmax -- the R equivalent of the Python
# script's log_pos() helper.
log_range <- log10(vmax / vmin)
log_pos <- function(value) log10(value / vmin) / log_range

# Same color stops as the other lightning scripts: white at vmin is a thin
# sliver, Tiffany blue holds as a plateau from ~0.0045 to 0.03, then
# yellow -> orange -> red up to 30.
lightning_colors <- c(
  '#ffffff', # white          (0.003, thin sliver)
  '#0ABAB5', # Tiffany blue starts (0.0045)
  '#0ABAB5', # Tiffany blue plateau ends (0.03)
  '#ffff00', # yellow         (0.3)
  '#ff8c00', # orange         (3)
  '#ff0000'  # red            (30)
)
lightning_positions <- log_pos(c(0.003, 0.0045, 0.03, 0.3, 3.0, 30.0))

make_custom_palette <- function(colors, positions, n) {
  rgb_mat  <- t(col2rgb(colors)) / 255
  eval_pos <- seq(0, 1, length.out = n)
  r <- approx(positions, rgb_mat[, 1], xout = eval_pos)$y
  g <- approx(positions, rgb_mat[, 2], xout = eval_pos)$y
  b <- approx(positions, rgb_mat[, 3], xout = eval_pos)$y
  rgb(r, g, b)
}

n_colors <- 60
lightning_palette <- make_custom_palette(lightning_colors, lightning_positions, n_colors)

# Log-spaced breaks (not linear) -- vmin * (vmax/vmin)^(i/n_colors), i.e.
# evenly spaced in log10-space, matching the Plate Carrée/Robinson scripts.
breaks <- vmin * (vmax / vmin) ^ (seq(0, n_colors) / n_colors)

# ── Natural Earth basemap layers ───────────────────────────────────────────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
# US only -- ne_states() with no country filter pulls every country's
# first-level admin boundaries (provinces, states, etc. worldwide), which
# produces dense internal hatching. We only want that level of detail for
# the US.
states    <- ne_states(country = "United States of America", returnclass = "sf")

# Reproject once, not per frame
coastline_p <- st_transform(coastline, target_crs)
borders_p   <- st_transform(borders, target_crs)
states_p    <- st_transform(states, target_crs)

# ── Load pre-computed climatology ──────────────────────────────────────────
# No rotate() needed: unlike ERA5, WGLC's longitude coordinate is already
# -180 to 180. No unit conversion either -- already annualized (strokes
# km-2 yr-1) by the climatology-build step.
r_stack <- rast(data_file, subds = "lightning_density_annual")
n_frames <- nlyr(r_stack)

nc <- nc_open(data_file)
frame_labels <- ncvar_get(nc, "frame_label")
nc_close(nc)

cat(sprintf("Loaded %d frames\n", n_frames))

# ── Helper: project a single lon/lat point into a target CRS ───────────────
proj_point <- function(lon, lat, target_crs) {
  pts <- vect(cbind(lon, lat), crs = "EPSG:4326")
  pts_proj <- project(pts, target_crs)
  crds(pts_proj)
}

# ── Helper: build a graticule line (meridian or parallel) in a target CRS ──
make_meridian <- function(lon, target_crs) {
  lats <- seq(-89.9, 89.9, length.out = 100)
  ln <- st_linestring(cbind(lon, lats))
  st_transform(st_sfc(ln, crs = "EPSG:4326"), target_crs)
}
make_parallel <- function(lat, target_crs) {
  lons <- seq(-180, 180, length.out = 200)
  ln <- st_linestring(cbind(lons, lat))
  st_transform(st_sfc(ln, crs = "EPSG:4326"), target_crs)
}

meridians <- lapply(c(-180, -90, 0, 90, 180), make_meridian, target_crs = target_crs)
parallels <- lapply(c(-90, -45, 0, 45, 90), make_parallel, target_crs = target_crs)

# ── Build clipping boundary (the true Foucaut diamond outline) ─────────────
# terra::project() on a raster fills the full rectangular bounding box and
# inverse-projects each cell back to a source lon/lat; outside the true
# projection domain that comes back wrapped/extrapolated rather than NA,
# producing a striped background outside the continents. This traces the
# up-the-antimeridian-and-back outline and masks the raster to it after
# projecting, matching the Python version's boundary handling.
lats_b    <- seq(-90, 90, length.out = 200)
left_pts  <- vect(cbind(-180, lats_b), crs = "EPSG:4326")
right_pts <- vect(cbind(180, rev(lats_b)), crs = "EPSG:4326")
boundary_pts    <- rbind(project(left_pts, target_crs), project(right_pts, target_crs))
boundary_coords <- crds(boundary_pts)
boundary_coords <- rbind(boundary_coords, boundary_coords[1, ])  # close the ring
boundary_poly   <- vect(boundary_coords, type = "polygons", crs = target_crs)

# ── Generate frames + MP4 ───────────────────────────────────────────────────
frame_paths <- character(n_frames)

for (i in seq_len(n_frames)) {
  
  field <- r_stack[[i]]
  
  # Log scale can't render exact zero (common over oceans/poles) and has
  # no "extend" mechanism the way ggplot2/matplotlib do -- clamp() both
  # ends BEFORE reprojecting, so bilinear interpolation during project()
  # never sees a raw exact-zero or extreme-outlier value.
  field_clamped <- clamp(field, lower = vmin, upper = vmax, values = TRUE)
  field_p <- mask(project(field_clamped, target_crs), boundary_poly)
  
  label <- frame_labels[i]
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(fpath, width = 1800, height = 980, res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1.3))
  
  # ── Map panel ──
  par(mar = c(2, 3, 7.5, 2), bty = "n")
  plot(field_p, col = lightning_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "")
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.6, col = "black")
  plot(st_geometry(borders_p), add = TRUE, lwd = 1.3, border = "black")
  plot(st_geometry(states_p), add = TRUE, lwd = 0.9, border = "black")
  
  for (m in meridians) plot(m, add = TRUE, col = "gray", lwd = 0.3)
  for (p in parallels) plot(p, add = TRUE, col = "gray", lwd = 0.3)
  
  # Pole labels -- Foucaut's poles are single points, not an edge to walk
  # longitudes along like Robinson/PlateCarree, so these are two fixed
  # annotations rather than a loop.
  xy_n <- proj_point(0, 90, target_crs)
  text(xy_n[1], xy_n[2], labels = "90\u00b0N", pos = 3, cex = 0.95,
       xpd = NA, offset = 0.3)
  xy_s <- proj_point(0, -90, target_crs)
  text(xy_s[1], xy_s[2], labels = "90\u00b0S", pos = 1, cex = 0.95,
       xpd = NA, offset = 0.3)
  
  # Latitude labels along the left edge only. No longitude labels --
  # Foucaut's pointed poles make bottom-edge longitude labels overlap into
  # unreadable text, same reasoning as the Python version.
  lat_labels <- list(c(-45, "45\u00b0S"), c(0, "0\u00b0"), c(45, "45\u00b0N"))
  for (ll in lat_labels) {
    lat_val <- as.numeric(ll[1]); lab <- ll[2]
    xy <- proj_point(-180, lat_val, target_crs)
    text(xy[1], xy[2], labels = lab, pos = 2, cex = 0.95, xpd = NA, offset = 0.4)
  }
  
  mtext("WWLLN | WGLC", side = 3, line = 6, adj = 0,
        cex = 1.0, font = 2)
  mtext(sprintf("%s; 2010\u20132025 Climatology", label), side = 3, line = 6,
        adj = 1, cex = 1.0, font = 2)
  
  # ── Colorbar panel (log-scale, edge-based image() call) ──────────────
  # x MUST be the uniform integer swatch index (0:n_colors), NOT the raw
  # log-spaced break VALUES -- image()'s x-axis is linear by default, so
  # passing raw log-spaced values as x would make each color cell's pixel
  # width proportional to the gap between consecutive break values, which
  # balloons hugely toward vmax (the 3->10 cell would render ~230x wider
  # on-screen than the 0.003->0.01 cell). Using the uniform index instead
  # gives equal pixel width per swatch, matching how the palette itself is
  # evenly spaced in log-FRACTION space -- the same "equal width per
  # decade" look as the PlateCarree/Robinson rect()-based colorbars, just
  # drawn via image() to match this script's established convention.
  par(mar = c(4.5, 6, 1.2, 6), bty = "n")
  cbar_mat <- matrix(seq_len(n_colors), ncol = 1)
  pad <- n_colors * 0.25
  
  image(x = 0:n_colors, y = 1, z = cbar_mat,
        col = lightning_palette, breaks = seq(0.5, n_colors + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(0 - pad, n_colors + pad))
  usr <- par("usr")
  rect(0, usr[3], n_colors, usr[4], border = "black", lwd = 1.2)
  
  # Tick positions use log_pos() * n_colors to convert each tick value into
  # its fractional position along the swatch index, then place the actual
  # tick_values as labels there.
  tick_pos <- log_pos(tick_values) * n_colors
  axis(1, at = tick_pos, labels = tick_values, cex.axis = 1.0)
  mtext(expression("Lightning stroke density (strokes km" ^ -2 * " yr" ^ -1 * ")"),
        side = 1, line = 3.0, cex = 1.0)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_path <- file.path(output_dir, "lightning_climatology_foucaut.mp4")
av_encode_video(frame_paths, output = video_path, framerate = 8, codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))