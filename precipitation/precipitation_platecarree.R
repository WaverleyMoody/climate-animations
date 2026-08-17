# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: precipitation_platecarree.R
# Description: Generates the precipitation climatology animation from ERA5 reanalysis data (1979-2000), rendered in the Plate Carrée projection. Weekly totals (cm water equivalent) computed server-side by the derived-era5-single-levels-daily-statistics CDS dataset, averaged into a 22-year climatological mean per 48-frame weekly cycle.
# Note: For the Robinson, Foucaut, and Nicolosi projections, see the other scripts in the precipitation scripts folder.

library(ncdf4)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────────
data_file  <- "/Volumes/CLIMATEDATA/precip/climatology/climatology_precip_48frame.nc"
output_dir <- "/Users/waverleymoody/Downloads/precipitation_animation/r_output"
frames_dir <- file.path(output_dir, "platecarree", "frames")
dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)

target_crs <- "+proj=longlat +datum=WGS84"

# ── Plot settings ─────────────────────────────────────────────────────────────
vmin <- 0
vmax <- 21
tick_step <- 3

# Same gradient as the Python Plate Carrée/Robinson/Foucaut/Nicolosi
# precipitation scripts: white at 0 cm blends into blue over just the
# first 0.3 cm, then continues through blue -> green -> yellow -> orange
# -> red -> purple up to 21 cm. Positions are fractions of vmax, mirroring
# the Python color_stops list exactly.
precip_colors <- c(
  '#ffffff', # white           (0 cm)
  '#d0e8f0', # pale blue       (0.3 cm)
  '#7ab8e0', # light blue      (2 cm)
  '#3060c0', # blue            (3 cm)
  '#00994d', # green           (6 cm)
  '#66cc00', # yellow-green    (8 cm)
  '#ffff00', # yellow          (9 cm)
  '#ff9900', # orange          (10 cm)
  '#ff3300', # red-orange      (11 cm)
  '#990000', # dark red        (12 cm)
  '#660066'  # purple          (21 cm)
)
precip_positions <- c(0, 0.3, 2, 3, 6, 8, 9, 10, 11, 12, 21) / vmax

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

n_colors      <- 60  # matches the Python script's 60 color steps
precip_palette <- make_custom_palette(precip_colors, precip_positions, n_colors)
breaks         <- seq(vmin, vmax, length.out = n_colors + 1)

# ── Natural Earth basemap layers ───────────────────────────────────────────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
states    <- ne_states(country = "United States of America", returnclass = "sf")

# Reproject once, not per frame (identity transform for Plate Carrée, but
# keeps this script structurally consistent with the Robinson/Foucaut
# versions where st_transform() actually does something).
coastline_p <- st_transform(coastline, target_crs)
borders_p   <- st_transform(borders, target_crs)
states_p    <- st_transform(states, target_crs)

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

# ── Load pre-computed climatology ──────────────────────────────────────────
r_stack <- rast(data_file, subds = "tp_weekly_total_cm")
r_stack <- rotate(r_stack)  # ERA5 is 0-360 -- fix before plotting
n_frames <- nlyr(r_stack)

nc <- nc_open(data_file)
frame_labels <- ncvar_get(nc, "frame_label")
nc_close(nc)

cat(sprintf("Loaded %d frames\n", n_frames))

# ── Generate PNG frames ──────────────────────────────────────────────────────
frame_paths <- character(n_frames)

for (i in seq_len(n_frames)) {
  
  field_p <- r_stack[[i]]  # already lon/lat after rotate(), no reprojection needed
  label <- frame_labels[i]
  
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(filename = fpath, width = 1800, height = 980, res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1))
  
  # ── Panel 1: map ──────────────────────────────────────────────────────
  par(mar = c(2, 3, 7.5, 2))
  
  plot(field_p, col = precip_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "", box = FALSE)
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.0, col = "black")
  plot(st_geometry(borders_p), add = TRUE, lwd = 0.8, border = "black")
  plot(st_geometry(states_p), add = TRUE, lwd = 0.5, border = "black")
  
  for (m in meridians) plot(m, add = TRUE, col = "gray", lwd = 0.3)
  for (p in parallels) plot(p, add = TRUE, col = "gray", lwd = 0.3)
  
  # ── Lat/lon axis labels (90°N/45°N/0°/45°S/90°S and 180°/90°W/0°/90°E/180°) ──
  lon_labels <- list(c(-180, "180\u00b0"), c(-90, "90\u00b0W"), c(0, "0\u00b0"),
                     c(90, "90\u00b0E"), c(180, "180\u00b0"))
  for (ll in lon_labels) {
    lon_val <- as.numeric(ll[1]); lab <- ll[2]
    text(lon_val, -95, labels = lab, pos = 1, cex = 0.75, xpd = NA)
  }
  
  lat_labels <- list(c(-90, "90\u00b0S"), c(-45, "45\u00b0S"), c(0, "0\u00b0"),
                     c(45, "45\u00b0N"), c(90, "90\u00b0N"))
  for (ll in lat_labels) {
    lat_val <- as.numeric(ll[1]); lab <- ll[2]
    text(-183, lat_val, labels = lab, pos = 2, cex = 0.75, xpd = NA)
  }
  
  # Title text, placed well above the 90N corner label.
  mtext("ERA-5 | Climate Reanalyzer", side = 3, line = 6, adj = 0,
        cex = 1.0, font = 2)
  mtext(sprintf("%s; 1979\u20132000 Weekly Total Mean", label), side = 3, line = 6,
        adj = 1, cex = 1.0, font = 2)
  
  # ── Panel 2: colorbar ─────────────────────────────────────────────────
  par(mar = c(3, 6, 0.5, 6), bty = "n")
  pad    <- (vmax - vmin) * 0.25
  cbar_z <- matrix(seq_len(n_colors), ncol = 1)
  image(x = breaks, y = c(0.5, 1.5), z = cbar_z,
        col = precip_palette, breaks = seq(0.5, n_colors + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(vmin - pad, vmax + pad))
  usr <- par("usr")
  rect(vmin, usr[3], vmax, usr[4], border = "black", lwd = 1.2)
  axis(1, at = seq(vmin, vmax, tick_step), labels = seq(vmin, vmax, tick_step),
       cex.axis = 1.0)
  mtext("Precipitation (cm)", side = 1, line = 2.0, cex = 1.0)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

# ── Assemble into MP4 ─────────────────────────────────────────────────────────
video_path <- file.path(output_dir, "precipitation_climatology.mp4")
av_encode_video(frame_paths, output = video_path, framerate = 4, codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))