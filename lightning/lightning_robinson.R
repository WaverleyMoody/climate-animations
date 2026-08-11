# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: lightning_Robinson.R
# Description: Generates the lightning climatology animation from WWLLN/WGLC data (2010-2025), rendered in the Robinson projection. 365 calendar-day frames (Feb 29 excluded), each the multi-year daily-climatological mean, annualized (x365) to match the reference product's strokes km-2 yr-1 units.
#
# Note: For the Plate Carrée, Foucaut, and Nicolosi projections, see the other scripts in the lightning scripts folder.


library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Volumes/CLIMATEDATA/lightning/climatology/climatology_lightning_365frame.nc"
output_dir <- "/Users/waverleymoody/Downloads/lightning_animation/r_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Projection to render ────────────────────────────────────────────────────
proj_name  <- "robinson"
target_crs <- "+proj=robin +datum=WGS84 +units=m +no_defs"

# ── Plot settings (log scale) ────────────────────────────────────────────
vmin <- 0.003
vmax <- 30
tick_values <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)

# log_pos(value) maps a data value to its fractional position (0-1) along
# the log10 range from vmin to vmax -- the R equivalent of the Python
# script's log_pos() helper, used identically for the palette color-stop
# positions and the colorbar tick positions.
log_range <- log10(vmax / vmin)
log_pos <- function(value) log10(value / vmin) / log_range

# Same color stops as the Python Plate Carrée/Robinson/Foucaut/Nicolosi
# lightning scripts: white at vmin is a thin sliver, Tiffany blue holds as
# a plateau from ~0.0045 to 0.03 (repeated at two anchor positions to force
# a flat span rather than an instant transition), then yellow -> orange ->
# red up to 30.
lightning_colors <- c(
  '#ffffff', # white          (0.003, thin sliver)
  '#0ABAB5', # Tiffany blue starts (0.0045)
  '#0ABAB5', # Tiffany blue plateau ends (0.03)
  '#ffff00', # yellow         (0.3)
  '#ff8c00', # orange         (3)
  '#ff0000'  # red            (30)
)
lightning_positions <- log_pos(c(0.003, 0.0045, 0.03, 0.3, 3.0, 30.0))

# Mirrors LinearSegmentedColormap.from_list with (position, color) tuples,
# since colorRampPalette() only supports evenly-spaced stops.
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

# Log-spaced breaks (not linear) -- the R equivalent of the Python script's
# np.logspace()-based log_levels/LogNorm. Each break is
# vmin * (vmax/vmin)^(i/n_colors), i.e. evenly spaced in log10-space.
breaks <- vmin * (vmax / vmin) ^ (seq(0, n_colors) / n_colors)

# ── Natural Earth basemap layers (cached locally after first download) ─────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
states    <- ne_states(country = "United States of America", returnclass = "sf")

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

# ── Generate frames + MP4 ───────────────────────────────────────────────────
cat(sprintf("=== Rendering projection: %s ===\n", proj_name))

frames_dir <- file.path(output_dir, proj_name, "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)
frame_paths <- character(n_frames)

# Reproject basemap + graticule layers once (not per frame)
coastline_p <- st_transform(coastline, target_crs)
borders_p   <- st_transform(borders, target_crs)
states_p    <- st_transform(states, target_crs)

grat_lons <- c(-180, -90, 0, 90, 180)
grat_lats <- c(-90, -45, 0, 45, 90)
meridians <- lapply(grat_lons, make_meridian, target_crs = target_crs)
parallels <- lapply(grat_lats, make_parallel, target_crs = target_crs)

for (i in seq_len(n_frames)) {
  field <- r_stack[[i]]
  
  # Log scale can't render exact zero (common over oceans/poles) and has
  # no "extend" mechanism the way ggplot2/matplotlib do -- clamp() both
  # ends so near-zero cells get the lowest color instead of NA, and cells
  # above vmax get the top color instead of NA (mirroring the Python
  # script's np.clip() + extend='max' behavior).
  field_clamped <- clamp(field, lower = vmin, upper = vmax, values = TRUE)
  field_p <- project(field_clamped, target_crs, method = "bilinear")
  
  label <- as.character(frame_labels[i])
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(fpath, width = 1800, height = 980, units = "px", res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1))
  
  # ── Panel 1: map ──────────────────────────────────────────────────────
  par(mar = c(2, 3, 7.5, 2))
  
  plot(field_p, col = lightning_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "", box = FALSE)
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.0, col = "black")
  plot(st_geometry(borders_p), add = TRUE, lwd = 0.8, border = "black")
  plot(st_geometry(states_p), add = TRUE, lwd = 0.5, border = "black")
  
  for (m in meridians) plot(m, add = TRUE, col = "gray", lwd = 0.3)
  for (p in parallels) plot(p, add = TRUE, col = "gray", lwd = 0.3)
  
  # Longitude labels along the bottom edge
  lon_labels <- list(c(-180, "180\u00b0"), c(-90, "90\u00b0W"), c(0, "0\u00b0"),
                     c(90, "90\u00b0E"), c(180, "180\u00b0"))
  for (ll in lon_labels) {
    lon_val <- as.numeric(ll[1]); lab <- ll[2]
    xy <- proj_point(lon_val, -90, target_crs)
    text(xy[1], xy[2], labels = lab, pos = 1, cex = 0.75, xpd = NA, offset = 0.3)
  }
  
  # Latitude labels along the left edge
  # Robinson gets extra offset at the poles so the label clears the curved edge.
  lat_labels <- list(c(-90, "90\u00b0S"), c(-45, "45\u00b0S"), c(0, "0\u00b0"),
                     c(45, "45\u00b0N"), c(90, "90\u00b0N"))
  for (ll in lat_labels) {
    lat_val <- as.numeric(ll[1]); lab <- ll[2]
    xy <- proj_point(-180, lat_val, target_crs)
    off <- if (abs(lat_val) == 90) 1.0 else 0.4
    text(xy[1], xy[2], labels = lab, pos = 2, cex = 0.75, xpd = NA, offset = off)
  }
  
  # Title text (top-left / top-right), placed well above the 90N corner label
  mtext("WWLLN | WGLC", side = 3, line = 6, adj = 0,
        cex = 0.9, font = 2)
  mtext(sprintf("%s; 2010\u20132025 Climatology", label), side = 3, line = 6,
        adj = 1, cex = 0.9, font = 2)
  
  # ── Panel 2: custom horizontal colorbar (log-scale tick placement) ────
  par(mar = c(3, 6, 0.5, 6))
  n_swatches <- length(lightning_palette)
  plot(NA, xlim = c(0, n_swatches), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rect(0:(n_swatches - 1), 0, 1:n_swatches, 1, col = lightning_palette, border = NA)
  box(col = "black", lwd = 1)
  
  # Tick positions use log_pos(), NOT the linear (val - vmin)/(vmax - vmin)
  # formula the temperature script uses -- the swatches are laid out evenly
  # in log-space (see make_custom_palette's eval_pos), so tick placement
  # has to follow the same log-space fraction to land correctly.
  tick_pos <- log_pos(tick_values) * n_swatches
  axis(1, at = tick_pos, labels = tick_values, cex.axis = 0.8, tck = -0.3,
       mgp = c(3, 0.5, 0))
  mtext(expression("Lightning stroke density (strokes km" ^ -2 * " yr" ^ -1 * ")"),
        side = 1, line = 2, cex = 0.9)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_name <- "lightning_climatology_robinson.mp4"
video_path <- file.path(output_dir, video_name)

av_encode_video(frame_paths, output = video_path, framerate = 8,
                codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))