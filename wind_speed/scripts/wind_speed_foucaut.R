# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: wind_speed_foucaut.R
# Description: Generates the 10m surface wind speed climatology animation from ERA5 reanalysis data (1979-2000), with mean wind vectors overlaid, rendered in the Foucaut projection.
#
# Note: For the Plate Carrée, Robinson, and Nicolosi projections, see the other scripts in the wind_speed scripts folder.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Users/waverleymoody/Downloads/climate_data_by_variable/wind_climatology.nc"
output_dir <- "/Users/waverleymoody/Downloads/wind_animation_R"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

proj_name  <- "foucaut"
target_crs <- "+proj=fouc +datum=WGS84 +units=m +no_defs"

frames_dir <- file.path(output_dir, proj_name, "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

# ── Plot settings (identical to animate_wind.py) ───────────────────────────
vmin <- 0
vmax <- 16

# Starting point from the Robinson/PlateCarree tuning session -- bumped up
# from 0.7 for longer arrow tails on Foucaut.
arrow_scale_deg <- 1.0
arrow_step      <- 35

# ── Custom colormap with uneven stop positions ──────────────────────────────
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

wind_color_stops <- c(
  '#ffffff', # white         (0 m/s)
  '#d0e8f0', # light blue    (3 m/s)
  '#469a46', # light green   (5.5 m/s)
  '#ffff00', # yellow        (7 m/s)
  '#ff6600', # orange        (8 m/s)
  '#ff3300', # darker orange (9.5 m/s)
  '#ff0000', # light red     (11 m/s)
  '#cc0000', # red           (13 m/s)
  '#990000', # dark red      (14.5 m/s)
  '#660000'  # darkest red   (16 m/s)
)
wind_positions <- c(0.00, 0.20, 0.35, 0.44, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

n_colors <- 35
wind_palette <- make_custom_palette(wind_color_stops, wind_positions, n_colors)
breaks <- seq(vmin, vmax, length.out = n_colors + 1)

# ── Natural Earth basemap layers ───────────────────────────────────────────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
states    <- ne_states(country = "United States of America", returnclass = "sf")

# Reproject once, not per frame
coastline_p <- st_transform(coastline, target_crs)
borders_p   <- st_transform(borders, target_crs)
states_p    <- st_transform(states, target_crs)

# ── Load pre-computed climatology ──────────────────────────────────────────
speed_stack <- rast(data_file, subds = "wind_speed")
u_stack     <- rast(data_file, subds = "u10")
v_stack     <- rast(data_file, subds = "v10")

# ERA5 is 0-360 -- fix all three consistently before reprojecting/subsampling
speed_stack <- rotate(speed_stack)
u_stack     <- rotate(u_stack)
v_stack     <- rotate(v_stack)
n_frames    <- nlyr(speed_stack)

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
# See 2m_temp_foucaut.R for the full explanation.
lats_b    <- seq(-90, 90, length.out = 200)
left_pts  <- vect(cbind(-180, lats_b), crs = "EPSG:4326")
right_pts <- vect(cbind(180, rev(lats_b)), crs = "EPSG:4326")
boundary_pts    <- rbind(project(left_pts, target_crs), project(right_pts, target_crs))
boundary_coords <- crds(boundary_pts)
boundary_coords <- rbind(boundary_coords, boundary_coords[1, ])  # close the ring
boundary_poly   <- vect(boundary_coords, type = "polygons", crs = target_crs)

# ── Arrow grid (lon/lat positions only -- frame-invariant, built once) ─────
col_idx  <- seq(1, ncol(u_stack), by = arrow_step)
row_idx  <- seq(1, nrow(u_stack), by = arrow_step)
lons_sub <- xFromCol(u_stack, col_idx)
lats_sub <- yFromRow(u_stack, row_idx)

idx_grid <- expand.grid(ri = seq_along(row_idx), ci = seq_along(col_idx))
grid_lon <- lons_sub[idx_grid$ci]
grid_lat <- lats_sub[idx_grid$ri]

# ── Generate frames + MP4 ───────────────────────────────────────────────────
frame_paths <- character(n_frames)

for (i in seq_len(n_frames)) {
  
  speed_field_p <- mask(project(speed_stack[[i]], target_crs), boundary_poly)
  
  # Per-frame u/v values at the fixed arrow grid positions
  u_mat   <- as.matrix(u_stack[[i]], wide = TRUE)
  v_mat   <- as.matrix(v_stack[[i]], wide = TRUE)
  grid_u  <- u_mat[cbind(row_idx[idx_grid$ri], col_idx[idx_grid$ci])]
  grid_v  <- v_mat[cbind(row_idx[idx_grid$ri], col_idx[idx_grid$ci])]
  
  end_lon <- grid_lon + grid_u * arrow_scale_deg
  end_lon <- ((end_lon + 180) %% 360) - 180   # dateline wrap
  end_lat <- pmin(pmax(grid_lat + grid_v * arrow_scale_deg, -90), 90)
  
  valid    <- is.finite(grid_u) & is.finite(grid_v)
  start_xy <- proj_point(grid_lon[valid], grid_lat[valid], target_crs)
  end_xy   <- proj_point(end_lon[valid], end_lat[valid], target_crs)
  
  # Outlier filter -- drops the handful of bad-projection points per frame
  # that produced stray horizontal-line artifacts on Robinson. This is the
  # fix that actually matters; arrow_scale_deg/arrow_step are cosmetic.
  map_ext       <- ext(speed_field_p)
  map_width     <- as.numeric(map_ext[2] - map_ext[1])
  max_arrow_len <- map_width * 0.03
  arrow_len     <- sqrt((end_xy[, 1] - start_xy[, 1])^2 + (end_xy[, 2] - start_xy[, 2])^2)
  keep          <- arrow_len < max_arrow_len
  
  label <- frame_labels[i]
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(fpath, width = 1800, height = 980, res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1.3))
  
  # ── Map panel ──
  par(mar = c(2, 3, 7.5, 2), bty = "n")
  plot(speed_field_p, col = wind_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "")
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.0, col = "black")
  plot(st_geometry(borders_p), add = TRUE, lwd = 0.8, border = "black")
  plot(st_geometry(states_p), add = TRUE, lwd = 0.5, border = "black")
  
  for (m in meridians) plot(m, add = TRUE, col = "gray", lwd = 0.3)
  for (p in parallels) plot(p, add = TRUE, col = "gray", lwd = 0.3)
  
  arrows(start_xy[keep, 1], start_xy[keep, 2], end_xy[keep, 1], end_xy[keep, 2],
         length = 0.03, angle = 20, code = 2, lwd = 0.6, col = "black")
  
  # Pole labels -- single points, not a longitude loop (same as temp/SLP)
  xy_n <- proj_point(0, 90, target_crs)
  text(xy_n[1], xy_n[2], labels = "90\u00b0N", pos = 3, cex = 0.95,
       xpd = NA, offset = 0.3)
  xy_s <- proj_point(0, -90, target_crs)
  text(xy_s[1], xy_s[2], labels = "90\u00b0S", pos = 1, cex = 0.95,
       xpd = NA, offset = 0.3)
  
  # Latitude labels along the left edge only. No longitude labels.
  lat_labels <- list(c(-45, "45\u00b0S"), c(0, "0\u00b0"), c(45, "45\u00b0N"))
  for (ll in lat_labels) {
    lat_val <- as.numeric(ll[1]); lab <- ll[2]
    xy <- proj_point(-180, lat_val, target_crs)
    text(xy[1], xy[2], labels = lab, pos = 2, cex = 0.95, xpd = NA, offset = 0.4)
  }
  
  mtext("ERA-5 | Climate Reanalyzer", side = 3, line = 6, adj = 0,
        cex = 1.0, font = 2)
  mtext(sprintf("%s; 1979\u20132000 Weekly Mean", label), side = 3, line = 6,
        adj = 1, cex = 1.0, font = 2)
  
  # ── Colorbar panel ──
  # No divider lines here -- animate_wind.py doesn't have them (unlike
  # animate_SLP.py), so none added to match.
  par(mar = c(4.5, 6, 1.2, 6), bty = "n")
  pad <- (vmax - vmin) * 0.25
  # x given as cell edges (breaks), not centers -- see SLP_foucaut.R for
  # why: centers let image() auto-extend the outer edges past vmin/vmax,
  # which poked color past the border rect.
  cbar_z <- matrix(seq_len(n_colors), ncol = 1)
  image(x = breaks, y = c(0.5, 1.5), z = cbar_z,
        col = wind_palette, breaks = seq(0.5, n_colors + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(vmin - pad, vmax + pad))
  usr <- par("usr")
  rect(vmin, usr[3], vmax, usr[4], border = "black", lwd = 1.2)
  axis(1, at = seq(vmin, vmax, 2), labels = seq(vmin, vmax, 2), cex.axis = 1.0)
  mtext("Wind Speed at 10 m (m/s)", side = 1, line = 3.0, cex = 1.0)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_path <- file.path(output_dir, "wind_climatology_foucaut.mp4")
av_encode_video(frame_paths, output = video_path, framerate = 4, codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))