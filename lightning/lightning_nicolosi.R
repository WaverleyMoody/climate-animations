# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: lightning_nicolosi.R
# Description: Generates the lightning climatology animation from WWLLN/WGLC 
# data (2010-2025), rendered as a double-hemisphere Nicolosi Globular projection. 
# 365 calendar-day frames (Feb 29 excluded), each the multi-year daily-climatological 
# mean, annualized (x365) to match the reference product's strokes km-2 yr-1 units.
#
# Note: For the Plate Carrée, Robinson, and Foucaut projections, see the other scripts 
# in the lightning scripts folder.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# Disable S2 spherical geometry for hemisphere clipping.
sf::sf_use_s2(FALSE)

# ── Paths ────────────────────────────────────────────────────────────────
data_file <- "/Volumes/CLIMATEDATA/lightning/climatology/climatology_lightning_365frame.nc"

output_dir <- "/Users/waverleymoody/Downloads/lightning_animation/r_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

frames_dir <- file.path(output_dir, "nicolosi", "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

# ── Projection settings ──────────────────────────────────────────────────
WEST_LON <- -90
EAST_LON <- 90

# Higher values create smoother output but increase processing time.
grid_res <- 600

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
  "#ffffff", # white          (0.003, thin sliver)
  "#0ABAB5", # Tiffany blue starts (0.0045)
  "#0ABAB5", # Tiffany blue plateau ends (0.03)
  "#ffff00", # yellow         (0.3)
  "#ff8c00", # orange         (3)
  "#ff0000"  # red            (30)
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

lightning_palette <- make_custom_palette(
  lightning_colors,
  lightning_positions,
  n_colors
)

# Log-spaced breaks (not linear) -- vmin * (vmax/vmin)^(i/n_colors), i.e.
# evenly spaced in log10-space. These feed terra::plot()'s data-value
# binning for the actual hemisphere rasters (unaffected by the
# linear-pixel-axis issue that applies specifically to the colorbar's
# image() call below).
breaks <- vmin * (vmax / vmin) ^ (seq(0, n_colors) / n_colors)

# ── Natural Earth basemap layers ─────────────────────────────────────────
cat("Downloading/loading basemap layers...\n")

coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
# US only -- ne_states() with no country filter pulls every country's
# first-level admin boundaries worldwide, which produces dense internal
# hatching. We only want that level of detail for the US.
states    <- ne_states(country = "United States of America", returnclass = "sf")

# ── Load pre-computed climatology ────────────────────────────────────────
# No rotate() needed: unlike ERA5, WGLC's longitude coordinate is already
# -180 to 180. No unit conversion either -- already annualized (strokes
# km-2 yr-1) by the climatology-build step.
r_stack <- rast(data_file, subds = "lightning_density_annual")
n_frames <- nlyr(r_stack)

lons <- xFromCol(r_stack, 1:ncol(r_stack))
lats <- yFromRow(r_stack, 1:nrow(r_stack))

nc <- nc_open(data_file)
frame_labels <- ncvar_get(nc, "frame_label")
nc_close(nc)

cat(sprintf("Loaded %d frames, grid %d x %d\n", n_frames, length(lats), length(lons)))

# ── Build one hemisphere ─────────────────────────────────────────────────
build_hemisphere <- function(lons, lats, central_lon, grid_res) {
  
  proj_str <- sprintf("+proj=nicol +lon_0=%f +R=6371000", central_lon)
  
  # Calculate angular longitude offset from the hemisphere center.
  lon_diff <- ((lons - central_lon + 180) %% 360) - 180
  
  # Keep only longitudes within 90 degrees of the center.
  mask_1d <- abs(lon_diff) <= 90
  candidate_idx <- which(mask_1d)
  
  # Sort longitudes by angular offset.
  ord <- order(lon_diff[candidate_idx])
  col_index <- candidate_idx[ord]
  
  lons_hemi <- lons[col_index]
  
  # Create the longitude/latitude grid. Latitude varies fastest, matching
  # R's column-major matrix flattening when as.vector() is used later.
  grid <- expand.grid(lat = lats, lon = lons_hemi)
  
  pts <- vect(cbind(grid$lon, grid$lat), crs = "EPSG:4326")
  XY <- crds(project(pts, proj_str))
  
  # Trace the outer boundary of the hemisphere.
  edge_lats <- seq(-90, 90, length.out = 400)
  
  b1 <- vect(cbind(rep(central_lon + 90, length(edge_lats)), edge_lats), crs = "EPSG:4326")
  b2 <- vect(cbind(rep(central_lon - 90, length(edge_lats)), edge_lats), crs = "EPSG:4326")
  
  b1_xy <- crds(project(b1, proj_str))
  b2_xy <- crds(project(b2, proj_str))
  
  boundary_xy <- rbind(b1_xy, b2_xy[nrow(b2_xy):1, ])
  boundary_poly <- vect(rbind(boundary_xy, boundary_xy[1, ]), type = "polygons", crs = proj_str)
  
  # Clip coastlines, borders, and states before projection.
  hemi_bbox <- st_as_sfc(
    st_bbox(
      c(xmin = central_lon - 90, xmax = central_lon + 90, ymin = -90, ymax = 90),
      crs = st_crs(4326)
    )
  )
  
  coastline_clip <- suppressWarnings(st_intersection(st_geometry(coastline), hemi_bbox))
  borders_clip   <- suppressWarnings(st_intersection(st_geometry(borders), hemi_bbox))
  states_clip    <- suppressWarnings(st_intersection(st_geometry(states), hemi_bbox))
  
  coastline_p <- st_transform(coastline_clip, proj_str)
  borders_p   <- st_transform(borders_clip, proj_str)
  states_p    <- st_transform(states_clip, proj_str)
  
  # Create the regular projected-space raster template.
  template <- rast(ext(boundary_poly), ncol = grid_res, nrow = grid_res, crs = proj_str)
  
  list(
    col_index = col_index,
    XY = XY,
    proj_str = proj_str,
    template = template,
    boundary_poly = boundary_poly,
    coastline_p = coastline_p,
    borders_p = borders_p,
    states_p = states_p
  )
}

cat("Precomputing Western Hemisphere geometry...\n")
west <- build_hemisphere(lons = lons, lats = lats, central_lon = WEST_LON, grid_res = grid_res)
cat("Precomputing Eastern Hemisphere geometry...\n")
east <- build_hemisphere(lons = lons, lats = lats, central_lon = EAST_LON, grid_res = grid_res)

# ── Rasterize one hemisphere for one animation frame ──────────────────────
render_hemisphere_raster <- function(frame_mat, hemi) {
  
  # Extract the data columns belonging to this hemisphere.
  data_hemi <- frame_mat[, hemi$col_index, drop = FALSE]
  
  # Flatten the matrix in column-major order.
  values_vec <- as.vector(data_hemi)
  
  # Log scale can't render exact zero (common over oceans/poles) and has
  # no "extend" mechanism the way ggplot2/matplotlib do -- clamp both ends
  # BEFORE rasterizing, so near-zero cells get the lowest color instead of
  # erroring, and cells above vmax get the top color instead of NA.
  values_vec <- pmin(pmax(values_vec, vmin), vmax)
  
  valid <- (
    is.finite(values_vec) &
      is.finite(hemi$XY[, 1]) &
      is.finite(hemi$XY[, 2])
  )
  
  pts_vals <- vect(hemi$XY[valid, , drop = FALSE], crs = hemi$proj_str)
  pts_vals$val <- values_vec[valid]
  
  # Rasterize the projected point values.
  r_hemi <- rasterize(pts_vals, hemi$template, field = "val", fun = "mean")
  
  # Fill empty cells caused by uneven projected point spacing.
  filled <- focal(r_hemi, w = 5, fun = "mean", na.rm = TRUE, na.policy = "only")
  r_hemi <- cover(r_hemi, filled)
  
  # Clip the result to the hemisphere boundary.
  mask(r_hemi, hemi$boundary_poly)
}

# ── Generate PNG frames ────────────────────────────────────────────────────
frame_paths <- character(n_frames)

for (i in seq_len(n_frames)) {
  
  frame_mat <- as.matrix(r_stack[[i]], wide = TRUE)
  
  west_r <- render_hemisphere_raster(frame_mat, west)
  east_r <- render_hemisphere_raster(frame_mat, east)
  
  label <- frame_labels[i]
  
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(filename = fpath, width = 1800, height = 980, res = 150)
  
  # Layout:
  #   Panel 1: full-width title row
  #   Panel 2: western hemisphere
  #   Panel 3: eastern hemisphere
  #   Panel 4: full-width colorbar
  layout(
    matrix(c(1, 1, 2, 3, 4, 4), nrow = 3, byrow = TRUE),
    heights = c(0.42, 5, 0.95)
  )
  
  # ── Full-width title panel ───────────────────────────────────────────
  par(mar = c(0, 1, 0, 1), bty = "n")
  plot.new()
  
  mtext("WWLLN | WGLC", side = 3, line = -1.8, adj = 0, cex = 1.0, font = 2)
  mtext(sprintf("%s; 2010\u20132025 Climatology", label), side = 3, line = -1.8,
        adj = 1, cex = 1.0, font = 2)
  
  # ── Western Hemisphere panel ─────────────────────────────────────────
  par(mar = c(1, 1, 1, 1), bty = "n")
  
  plot(west_r, col = lightning_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "", asp = 1)
  
  plot(west$coastline_p, add = TRUE, lwd = 1.6, col = "black")
  plot(west$borders_p, add = TRUE, lwd = 1.3, border = "black")
  plot(west$states_p, add = TRUE, lwd = 0.9, border = "black")
  plot(west$boundary_poly, add = TRUE, border = "black", lwd = 1, col = NA)
  
  # ── Eastern Hemisphere panel ─────────────────────────────────────────
  par(mar = c(1, 1, 1, 1), bty = "n")
  
  plot(east_r, col = lightning_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "", asp = 1)
  
  plot(east$coastline_p, add = TRUE, lwd = 1.6, col = "black")
  plot(east$borders_p, add = TRUE, lwd = 1.3, border = "black")
  plot(east$states_p, add = TRUE, lwd = 0.9, border = "black")
  plot(east$boundary_poly, add = TRUE, border = "black", lwd = 1, col = NA)
  
  # ── Full-width colorbar panel (log-scale, uniform swatch index) ──────
  par(mar = c(4.0, 6, 0.8, 6), bty = "n", pty = "m")
  
  cbar_z <- matrix(seq_len(n_colors), ncol = 1)
  pad <- n_colors * 0.25
  
  # x = 0:n_colors (uniform integer index), NOT raw log-spaced break
  # values -- see header note.
  image(
    x = 0:n_colors, y = c(0.7, 1.3), z = cbar_z,
    col = lightning_palette,
    breaks = seq(0.5, n_colors + 0.5, by = 1),
    axes = FALSE, xlab = "", ylab = "",
    xlim = c(0 - pad, n_colors + pad),
    ylim = c(0.5, 1.5)
  )
  
  rect(0, 0.7, n_colors, 1.3, border = "black", lwd = 1.2)
  
  tick_pos <- log_pos(tick_values) * n_colors
  axis(side = 1, at = tick_pos, labels = tick_values, cex.axis = 1.0)
  
  mtext(expression("Lightning stroke density (strokes km" ^ -2 * " yr" ^ -1 * ")"),
        side = 1, line = 2.7, cex = 1.0)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

# ── Encode frames as MP4 ────────────────────────────────────────────────────
video_path <- file.path(output_dir, "lightning_climatology_nicolosi.mp4")

av_encode_video(input = frame_paths, output = video_path, framerate = 8, codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))