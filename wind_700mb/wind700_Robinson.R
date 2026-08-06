# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: wind700_Robinson.R
# Description: Generates the 700mb wind speed climatology animation from
# ERA5 reanalysis data (1979-2000), with mean wind vectors overlaid,
# rendered in the Robinson projection.
#
# Note: For the Plate Carrée, Foucaut, and Nicolosi projections, see the
# other scripts in the wind_700mb scripts folder.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Volumes/CLIMATEDATA/wind700/climatology/climatology_wind700_48frame.nc"
output_dir <- "/Users/waverleymoody/Downloads/wind700_animation_R"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Projection to render ────────────────────────────────────────────────────
proj_name  <- "robinson"
target_crs <- "+proj=robin +datum=WGS84 +units=m +no_defs"

# ── Plot settings ────────────────────────────────────────────────────────
vmin <- 0
vmax <- 30

# arrow_scale_deg controls how many degrees of lon/lat displacement each
# 1 m/s of wind produces. 700mb winds are faster than 10m surface winds,
# so this is tuned smaller than the 10m script's 0.7 to keep arrows a
# similar visual length despite higher wind speeds.
arrow_scale_deg <- 0.5
arrow_step      <- 35

# ── Custom colormap with uneven stop positions ──────────────────────────────
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
  '#d0e8f0', # light blue    (6 m/s)
  '#469a46', # light green   (10.5 m/s)
  '#ffff00', # yellow        (13 m/s)
  '#ff6600', # orange        (15 m/s)
  '#ff3300', # darker orange (18 m/s)
  '#ff0000', # light red     (21 m/s)
  '#cc0000', # red           (24 m/s)
  '#990000', # dark red      (27 m/s)
  '#660000'  # darkest red   (30 m/s)
)
wind_positions <- c(0.00, 0.20, 0.35, 0.44, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

n_colors <- 35
wind_palette <- make_custom_palette(wind_color_stops, wind_positions, n_colors)
breaks <- seq(vmin, vmax, length.out = n_colors + 1)

# ── Natural Earth basemap layers ─────────────────────────────────────────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
states    <- ne_states(country = "United States of America", returnclass = "sf")

# ── Load pre-computed climatology ──────────────────────────────────────────
speed_stack <- rast(data_file, subds = "wind_speed_700mb")
u_stack     <- rast(data_file, subds = "u700")
v_stack     <- rast(data_file, subds = "v700")

# ERA5 uses 0-360 longitude; convert all three to -180/180 consistently.
speed_stack <- rotate(speed_stack)
u_stack     <- rotate(u_stack)
v_stack     <- rotate(v_stack)

n_frames <- nlyr(speed_stack)

# No frame_label variable -- generate week labels from frame index.
# 48 frames = 12 months x 4 weeks (starting on 1st, 8th, 15th, 22nd).
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
week_starts <- c(1, 8, 15, 22)
frame_labels <- as.vector(outer(week_starts, month_names,
                                function(d, m) paste(m, d)))[
                                  order(rep(seq_along(month_names), each = length(week_starts)))]

cat(sprintf("Loaded %d frames\n", n_frames))

# ── Precompute the subsampled arrow grid (same for every frame) ────────────
col_idx  <- seq(1, ncol(u_stack), by = arrow_step)
row_idx  <- seq(1, nrow(u_stack), by = arrow_step)
lons_sub <- xFromCol(u_stack, col_idx)
lats_sub <- yFromRow(u_stack, row_idx)

idx_grid <- expand.grid(ri = seq_along(row_idx), ci = seq_along(col_idx))
grid_lon <- lons_sub[idx_grid$ci]
grid_lat <- lats_sub[idx_grid$ri]

# ── Helper: project one or many lon/lat points into a target CRS ───────────
proj_point <- function(lon, lat, target_crs) {
  pts <- vect(cbind(lon, lat), crs = "EPSG:4326")
  pts_proj <- project(pts, target_crs)
  crds(pts_proj)
}

# ── Helper: build a graticule line (meridian or parallel) in a target CRS ──
# Built point-by-point via terra (same method as proj_point), not via
# sf::st_transform on a LINESTRING spanning the antimeridian -- that approach
# produced degenerate straight-line artifacts in Robinson.
make_meridian <- function(lon, target_crs) {
  lats <- seq(-89.9, 89.9, length.out = 100)
  proj_point(rep(lon, length(lats)), lats, target_crs)
}

make_parallel <- function(lat, target_crs) {
  lons <- seq(-179.9, 179.9, length.out = 200)
  proj_point(lons, rep(lat, length(lons)), target_crs)
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
  speed_field   <- speed_stack[[i]]
  speed_field_p <- project(speed_field, target_crs, method = "bilinear")
  
  label <- frame_labels[i]
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  # Pull u/v at the subsampled grid for this frame
  u_mat  <- as.matrix(u_stack[[i]], wide = TRUE)
  v_mat  <- as.matrix(v_stack[[i]], wide = TRUE)
  grid_u <- u_mat[cbind(row_idx[idx_grid$ri], col_idx[idx_grid$ci])]
  grid_v <- v_mat[cbind(row_idx[idx_grid$ri], col_idx[idx_grid$ci])]
  
  png(fpath, width = 1800, height = 980, units = "px", res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1))
  
  # ── Panel 1: map ──────────────────────────────────────────────────────
  par(mar = c(2, 3, 7.5, 2))
  
  plot(speed_field_p, col = wind_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "", box = FALSE)
  
  # Wind vector arrows: approximate the vector as a small lon/lat
  # displacement, then project both the start and end point into the
  # map CRS before drawing.
  end_lon <- grid_lon + grid_u * arrow_scale_deg
  end_lon <- ((end_lon + 180) %% 360) - 180
  end_lat <- pmin(pmax(grid_lat + grid_v * arrow_scale_deg, -90), 90)
  
  valid <- is.finite(grid_u) & is.finite(grid_v)
  start_xy <- proj_point(grid_lon[valid], grid_lat[valid], target_crs)
  end_xy   <- proj_point(end_lon[valid], end_lat[valid], target_crs)
  
  # Drop rare outlier arrows whose projected length is implausibly long
  # (edge-case points near the dateline/poles that produce long straight-
  # line artifacts across the map, especially in Robinson).
  map_ext       <- ext(speed_field_p)
  map_width     <- as.numeric(map_ext[2] - map_ext[1])
  max_arrow_len <- map_width * 0.03
  arrow_len     <- sqrt((end_xy[, 1] - start_xy[, 1])^2 +
                          (end_xy[, 2] - start_xy[, 2])^2)
  keep <- arrow_len < max_arrow_len
  
  arrows(start_xy[keep, 1], start_xy[keep, 2],
         end_xy[keep, 1],   end_xy[keep, 2],
         length = 0.03, angle = 20, code = 2, lwd = 0.6, col = "black")
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.0, col = "black")
  plot(st_geometry(borders_p),   add = TRUE, lwd = 0.8, border = "black")
  plot(st_geometry(states_p),    add = TRUE, lwd = 0.5, border = "black")
  
  for (m in meridians) lines(m[, 1], m[, 2], col = "gray", lwd = 0.3)
  for (p in parallels) lines(p[, 1], p[, 2], col = "gray", lwd = 0.3)
  
  # Longitude labels along the bottom edge
  lon_labels <- list(c(-180, "180\u00b0"), c(-90, "90\u00b0W"), c(0, "0\u00b0"),
                     c(90, "90\u00b0E"), c(180, "180\u00b0"))
  for (ll in lon_labels) {
    lon_val <- as.numeric(ll[1]); lab <- ll[2]
    xy <- proj_point(lon_val, -90, target_crs)
    text(xy[1], xy[2], labels = lab, pos = 1, cex = 0.75, xpd = NA, offset = 0.3)
  }
  
  # Latitude labels along the left edge.
  # Robinson gets extra offset at the poles so the label clears the curved edge.
  lat_labels <- list(c(-90, "90\u00b0S"), c(-45, "45\u00b0S"), c(0, "0\u00b0"),
                     c(45, "45\u00b0N"), c(90, "90\u00b0N"))
  for (ll in lat_labels) {
    lat_val <- as.numeric(ll[1]); lab <- ll[2]
    xy <- proj_point(-180, lat_val, target_crs)
    off <- if (abs(lat_val) == 90) 1.0 else 0.4
    text(xy[1], xy[2], labels = lab, pos = 2, cex = 0.75, xpd = NA, offset = off)
  }
  
  mtext("ERA-5 | Climate Reanalyzer", side = 3, line = 6, adj = 0,
        cex = 0.9, font = 2)
  mtext(sprintf("%s; 1979\u20132000 Weekly Mean", label), side = 3, line = 6,
        adj = 1, cex = 0.9, font = 2)
  
  # ── Panel 2: custom horizontal colorbar ──────────────────────────────
  par(mar = c(3, 6, 0.5, 6))
  n_swatches <- length(wind_palette)
  plot(NA, xlim = c(0, n_swatches), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rect(0:(n_swatches - 1), 0, 1:n_swatches, 1, col = wind_palette, border = NA)
  box(col = "black", lwd = 1)
  
  tick_vals <- seq(vmin, vmax, 5)
  tick_pos  <- (tick_vals - vmin) / (vmax - vmin) * n_swatches
  axis(1, at = tick_pos, labels = tick_vals, cex.axis = 0.8, tck = -0.3,
       mgp = c(3, 0.5, 0))
  mtext("Wind Speed at 700 mb (m/s)", side = 1, line = 2, cex = 0.9)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_name <- "wind700_climatology_robinson.mp4"
video_path <- file.path(output_dir, video_name)

av_encode_video(frame_paths, output = video_path, framerate = 4,
                codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))