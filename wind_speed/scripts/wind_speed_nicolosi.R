# SDSU Climate Informatics
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# San Diego State University
#
# R translation of animate_nicolosi_wind.py
# Full 48-frame double-hemisphere Nicolosi Globular animation for 10m
# wind speed, with mean wind vectors overlaid.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# Disable S2 spherical geometry for hemisphere clipping.
sf::sf_use_s2(FALSE)

# ── Paths ────────────────────────────────────────────────────────────────

data_file <- paste0(
  "/Users/waverleymoody/Downloads/",
  "climate_data_by_variable/wind_climatology.nc"
)

output_dir <- paste0(
  "/Users/waverleymoody/Downloads/",
  "nicolosi_wind_R"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

frames_dir <- file.path(
  output_dir,
  "frames"
)

dir.create(
  frames_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ── Projection settings ──────────────────────────────────────────────────

WEST_LON <- -90
EAST_LON <- 90

# Higher values create smoother output but increase processing time.
grid_res <- 600

# ── Quiver settings ──────────────────────────────────────────────────────

quiver_step <- 28

max_quiver_lat <- 85

# arrow_scale is in meters per (m/s)
arrow_scale <- 75000

# Secondary safety net -- caps arrow length as a fraction of the
# hemisphere's projected width
max_arrow_len_frac <- 0.35

arrow_lwd <- 0.6
arrow_head_length <- 0.025

# ── Plot settings ────────────────────────────────────────────────────────

# Same fixed color scale and custom colormap as wind_speed_platecarree.R /
# wind_speed_robinson.R, so this animation matches the PlateCarree/
# Robinson ones on the site.
vmin <- 0
vmax <- 16

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

wind_palette <- make_custom_palette(
  wind_color_stops,
  wind_positions,
  n_colors
)

breaks <- seq(
  vmin,
  vmax,
  length.out = n_colors + 1
)

# ── Natural Earth basemap layers ─────────────────────────────────────────

cat("Downloading/loading basemap layers...\n")

coastline <- ne_coastline(
  scale = 50,
  returnclass = "sf"
)

borders <- ne_countries(
  scale = 50,
  returnclass = "sf"
)

# ── Load pre-computed climatology ────────────────────────────────────────

speed_stack <- rast(data_file, subds = "wind_speed")
u_stack     <- rast(data_file, subds = "u10")
v_stack     <- rast(data_file, subds = "v10")

# ERA5 uses 0-360 longitude; convert all three to -180/180 consistently
speed_stack <- rotate(speed_stack)
u_stack     <- rotate(u_stack)
v_stack     <- rotate(v_stack)

n_frames <- nlyr(speed_stack)

lons <- xFromCol(
  speed_stack,
  1:ncol(speed_stack)
)

lats <- yFromRow(
  speed_stack,
  1:nrow(speed_stack)
)

nc <- nc_open(
  data_file
)

frame_labels <- ncvar_get(
  nc,
  "frame_label"
)

nc_close(nc)

cat(
  sprintf(
    "Loaded %d frames, grid %d x %d\n",
    n_frames,
    length(lats),
    length(lons)
  )
)

# ── Helper: project one or many lon/lat points into a target CRS ───────────
proj_point <- function(lon, lat, target_crs) {
  pts <- vect(cbind(lon, lat), crs = "EPSG:4326")
  pts_proj <- project(pts, target_crs)
  crds(pts_proj)
}

# ── Build one hemisphere ─────────────────────────────────────────────────

build_hemisphere <- function(
    lons,
    lats,
    central_lon,
    grid_res
) {
  
  proj_str <- sprintf(
    "+proj=nicol +lon_0=%f +R=6371000",
    central_lon
  )
  
  # Calculate angular longitude offset from the hemisphere center.
  lon_diff <- (
    (lons - central_lon + 180) %% 360
  ) - 180
  
  # Keep only longitudes within 90 degrees of the center.
  mask_1d <- abs(lon_diff) <= 90
  
  candidate_idx <- which(
    mask_1d
  )
  
  # Sort longitudes by angular offset.
  ord <- order(
    lon_diff[candidate_idx]
  )
  
  col_index <- candidate_idx[ord]
  
  lons_hemi <- lons[col_index]
  
  # Create the longitude/latitude grid.
  #
  # Latitude varies fastest, matching R's column-major matrix flattening
  # when as.vector() is used later.
  grid <- expand.grid(
    lat = lats,
    lon = lons_hemi
  )
  
  pts <- vect(
    cbind(
      grid$lon,
      grid$lat
    ),
    crs = "EPSG:4326"
  )
  
  XY <- crds(
    project(
      pts,
      proj_str
    )
  )
  
  # Trace the outer boundary of the hemisphere.
  edge_lats <- seq(
    -90,
    90,
    length.out = 400
  )
  
  b1 <- vect(
    cbind(
      rep(
        central_lon + 90,
        length(edge_lats)
      ),
      edge_lats
    ),
    crs = "EPSG:4326"
  )
  
  b2 <- vect(
    cbind(
      rep(
        central_lon - 90,
        length(edge_lats)
      ),
      edge_lats
    ),
    crs = "EPSG:4326"
  )
  
  b1_xy <- crds(
    project(
      b1,
      proj_str
    )
  )
  
  b2_xy <- crds(
    project(
      b2,
      proj_str
    )
  )
  
  boundary_xy <- rbind(
    b1_xy,
    b2_xy[nrow(b2_xy):1, ]
  )
  
  boundary_poly <- vect(
    rbind(
      boundary_xy,
      boundary_xy[1, ]
    ),
    type = "polygons",
    crs = proj_str
  )
  
  # Clip coastlines and borders before projection.
  hemi_bbox <- st_as_sfc(
    st_bbox(
      c(
        xmin = central_lon - 90,
        xmax = central_lon + 90,
        ymin = -90,
        ymax = 90
      ),
      crs = st_crs(4326)
    )
  )
  
  coastline_clip <- suppressWarnings(
    st_intersection(
      st_geometry(coastline),
      hemi_bbox
    )
  )
  
  borders_clip <- suppressWarnings(
    st_intersection(
      st_geometry(borders),
      hemi_bbox
    )
  )
  
  coastline_p <- st_transform(
    coastline_clip,
    proj_str
  )
  
  borders_p <- st_transform(
    borders_clip,
    proj_str
  )
  
  # Create the regular projected-space raster template.
  template <- rast(
    ext(boundary_poly),
    ncol = grid_res,
    nrow = grid_res,
    crs = proj_str
  )
  
  # Width of the hemisphere in projected meters -- used by the outlier
  # arrow-length filter in draw_hemisphere_quiver().
  boundary_ext <- ext(boundary_poly)
  map_width <- as.numeric(boundary_ext[2] - boundary_ext[1])
  
  # ── Quiver geometry: subsample the (already hemisphere-restricted,
  # correctly-ordered) row/col indices, project those points, and
  # precompute the local rotation Jacobian at each -- none of this
  # depends on frame data, only on the grid and the projection.
  quiver_rows_all <- seq(1, length(lats), by = quiver_step)
  
  # Drop rows too close to the pole -- see max_quiver_lat comment above.
  # This is a uniform latitude cutoff, so it's safe to filter entire rows
  # here rather than needing to touch quiver_cols.
  quiver_rows <- quiver_rows_all[abs(lats[quiver_rows_all]) <= max_quiver_lat]
  
  quiver_cols <- col_index[seq(1, length(col_index), by = quiver_step)]
  
  qgrid <- expand.grid(
    lat = lats[quiver_rows],
    lon = lons[quiver_cols]
  )
  
  quiver_XY <- proj_point(qgrid$lon, qgrid$lat, proj_str)
  
  # Use centered finite differences to estimate how local eastward and
  # northward vectors rotate under the Nicolosi projection. Centered
  # differences are more stable than one-sided differences near the curved
  # hemisphere edge.
  R_earth <- 6371000
  eps <- 0.01
  
  xy_east <- proj_point(
    qgrid$lon + eps,
    qgrid$lat,
    proj_str
  )
  
  xy_west <- proj_point(
    qgrid$lon - eps,
    qgrid$lat,
    proj_str
  )
  
  xy_north <- proj_point(
    qgrid$lon,
    pmin(qgrid$lat + eps, 89.999),
    proj_str
  )
  
  xy_south <- proj_point(
    qgrid$lon,
    pmax(qgrid$lat - eps, -89.999),
    proj_str
  )
  
  # Physical distances represented by the centered perturbations.
  east_dist <- (
    2 * eps * pi / 180 *
      R_earth *
      cos(qgrid$lat * pi / 180)
  )
  
  north_dist <- 2 * eps * pi / 180 * R_earth
  
  east_dist[abs(east_dist) < 1e-6] <- NA_real_
  
  # Local projection Jacobian:
  #
  #   [dx/de  dx/dn]
  #   [dy/de  dy/dn]
  #
  # where e and n are physical eastward and northward distances.
  quiver_dxde <- (
    xy_east[, 1] - xy_west[, 1]
  ) / east_dist
  
  quiver_dyde <- (
    xy_east[, 2] - xy_west[, 2]
  ) / east_dist
  
  quiver_dxdn <- (
    xy_north[, 1] - xy_south[, 1]
  ) / north_dist
  
  quiver_dydn <- (
    xy_north[, 2] - xy_south[, 2]
  ) / north_dist
  
  list(
    col_index = col_index,
    XY = XY,
    proj_str = proj_str,
    template = template,
    boundary_poly = boundary_poly,
    map_width = map_width,
    coastline_p = coastline_p,
    borders_p = borders_p,
    quiver_rows = quiver_rows,
    quiver_cols = quiver_cols,
    quiver_XY = quiver_XY,
    quiver_dxde = quiver_dxde,
    quiver_dyde = quiver_dyde,
    quiver_dxdn = quiver_dxdn,
    quiver_dydn = quiver_dydn
  )
}

cat("Precomputing Western Hemisphere geometry...\n")

west <- build_hemisphere(
  lons = lons,
  lats = lats,
  central_lon = WEST_LON,
  grid_res = grid_res
)

cat("Precomputing Eastern Hemisphere geometry...\n")

east <- build_hemisphere(
  lons = lons,
  lats = lats,
  central_lon = EAST_LON,
  grid_res = grid_res
)

# ── Rasterize one hemisphere for one animation frame ─────────────────────

render_hemisphere_raster <- function(
    frame_mat,
    hemi
) {
  
  # Extract the data columns belonging to this hemisphere.
  data_hemi <- frame_mat[
    ,
    hemi$col_index,
    drop = FALSE
  ]
  
  # Flatten the matrix in column-major order.
  values_vec <- as.vector(
    data_hemi
  )
  
  valid <- (
    is.finite(values_vec) &
      is.finite(hemi$XY[, 1]) &
      is.finite(hemi$XY[, 2])
  )
  
  pts_vals <- vect(
    hemi$XY[
      valid,
      ,
      drop = FALSE
    ],
    crs = hemi$proj_str
  )
  
  pts_vals$val <- values_vec[valid]
  
  # Rasterize the projected point values.
  r_hemi <- rasterize(
    pts_vals,
    hemi$template,
    field = "val",
    fun = "mean"
  )
  
  # Fill empty cells caused by uneven projected point spacing.
  filled <- focal(
    r_hemi,
    w = 5,
    fun = "mean",
    na.rm = TRUE,
    na.policy = "only"
  )
  
  r_hemi <- cover(
    r_hemi,
    filled
  )
  
  # Clip the result to the hemisphere boundary.
  mask(
    r_hemi,
    hemi$boundary_poly
  )
}

# ── Draw rotated wind vectors for one hemisphere on the current device ──

draw_hemisphere_quiver <- function(
    u_mat,
    v_mat,
    hemi
) {
  
  # Extract u and v in the same latitude-fastest ordering used when the
  # quiver coordinates and projection Jacobian were constructed.
  u_sub <- u_mat[
    hemi$quiver_rows,
    hemi$quiver_cols,
    drop = FALSE
  ]
  
  v_sub <- v_mat[
    hemi$quiver_rows,
    hemi$quiver_cols,
    drop = FALSE
  ]
  
  grid_u <- as.vector(u_sub)
  grid_v <- as.vector(v_sub)
  
  # ERA5 u/v components point toward the direction of air motion:
  #   u > 0 means eastward motion
  #   v > 0 means northward motion
  #
  # Transform the local east/north vector into projected x/y coordinates.
  dx_raw <- (
    hemi$quiver_dxde * grid_u +
      hemi$quiver_dxdn * grid_v
  )
  
  dy_raw <- (
    hemi$quiver_dyde * grid_u +
      hemi$quiver_dydn * grid_v
  )
  
  wind_speed <- sqrt(
    grid_u^2 + grid_v^2
  )
  
  projected_norm <- sqrt(
    dx_raw^2 + dy_raw^2
  )
  
  valid <- (
    is.finite(grid_u) &
      is.finite(grid_v) &
      is.finite(wind_speed) &
      is.finite(dx_raw) &
      is.finite(dy_raw) &
      is.finite(projected_norm) &
      projected_norm > 1e-10 &
      wind_speed > 0.05
  )
  
  start_xy <- hemi$quiver_XY[
    valid,
    ,
    drop = FALSE
  ]
  
  # Keep only the projected direction from the Jacobian, then restore the
  # original wind-speed magnitude. This prevents projection scale distortion
  # from producing misleading arrow lengths and apparent directions.
  direction_x <- dx_raw[valid] / projected_norm[valid]
  direction_y <- dy_raw[valid] / projected_norm[valid]
  
  display_length <- wind_speed[valid] * arrow_scale
  
  end_xy <- cbind(
    start_xy[, 1] + direction_x * display_length,
    start_xy[, 2] + direction_y * display_length
  )
  
  arrow_len <- sqrt(
    (end_xy[, 1] - start_xy[, 1])^2 +
      (end_xy[, 2] - start_xy[, 2])^2
  )
  
  max_arrow_len <- hemi$map_width * max_arrow_len_frac
  
  keep <- (
    is.finite(arrow_len) &
      arrow_len > 0 &
      arrow_len < max_arrow_len
  )
  
  arrows(
    x0 = start_xy[keep, 1],
    y0 = start_xy[keep, 2],
    x1 = end_xy[keep, 1],
    y1 = end_xy[keep, 2],
    length = arrow_head_length,
    angle = 20,
    code = 2,
    lwd = arrow_lwd,
    col = "black",
    xpd = FALSE
  )
}

# ── Generate PNG frames ──────────────────────────────────────────────────

frame_paths <- character(
  n_frames
)

for (i in seq_len(n_frames)) {
  
  speed_mat <- as.matrix(
    speed_stack[[i]],
    wide = TRUE
  )
  
  u_mat <- as.matrix(
    u_stack[[i]],
    wide = TRUE
  )
  
  v_mat <- as.matrix(
    v_stack[[i]],
    wide = TRUE
  )
  
  west_r <- render_hemisphere_raster(
    speed_mat,
    west
  )
  
  east_r <- render_hemisphere_raster(
    speed_mat,
    east
  )
  
  label <- frame_labels[i]
  
  fpath <- file.path(
    frames_dir,
    sprintf(
      "frame_%03d.png",
      i - 1
    )
  )
  
  png(
    filename = fpath,
    width = 1800,
    height = 980,
    res = 150
  )
  
  # Layout:
  #
  #   Panel 1: full-width title row
  #   Panel 2: western hemisphere
  #   Panel 3: eastern hemisphere
  #   Panel 4: full-width colorbar
  #
  # The first value controls title-row height.
  # The third value controls colorbar-row height.
  layout(
    matrix(
      c(
        1, 1,
        2, 3,
        4, 4
      ),
      nrow = 3,
      byrow = TRUE
    ),
    heights = c(
      0.42,
      5,
      0.95
    )
  )
  
  # ── Full-width title panel ─────────────────────────────────────────────
  
  par(
    mar = c(0, 1, 0, 1),
    bty = "n"
  )
  
  plot.new()
  
  # Left-side title.
  mtext(
    "ERA-5 | Climate Reanalyzer",
    side = 3,
    line = -1.8,
    adj = 0,
    cex = 1.0,
    font = 2
  )
  
  # Right-side title.
  mtext(
    sprintf(
      "%s; 1979\u20132000 Weekly Mean",
      label
    ),
    side = 3,
    line = -1.8,
    adj = 1,
    cex = 1.0,
    font = 2
  )
  
  # ── Western Hemisphere panel ───────────────────────────────────────────
  
  par(
    mar = c(1, 1, 1, 1),
    bty = "n"
  )
  
  plot(
    west_r,
    col = wind_palette,
    breaks = breaks,
    range = c(vmin, vmax),
    axes = FALSE,
    legend = FALSE,
    main = "",
    asp = 1
  )
  
  draw_hemisphere_quiver(
    u_mat,
    v_mat,
    west
  )
  
  plot(
    west$coastline_p,
    add = TRUE,
    lwd = 1.0,
    col = "black"
  )
  
  plot(
    west$borders_p,
    add = TRUE,
    lwd = 0.8,
    border = "black"
  )
  
  plot(
    west$boundary_poly,
    add = TRUE,
    border = "black",
    lwd = 1,
    col = NA
  )
  
  # ── Eastern Hemisphere panel ───────────────────────────────────────────
  
  par(
    mar = c(1, 1, 1, 1),
    bty = "n"
  )
  
  plot(
    east_r,
    col = wind_palette,
    breaks = breaks,
    range = c(vmin, vmax),
    axes = FALSE,
    legend = FALSE,
    main = "",
    asp = 1
  )
  
  draw_hemisphere_quiver(
    u_mat,
    v_mat,
    east
  )
  
  plot(
    east$coastline_p,
    add = TRUE,
    lwd = 1.0,
    col = "black"
  )
  
  plot(
    east$borders_p,
    add = TRUE,
    lwd = 0.8,
    border = "black"
  )
  
  plot(
    east$boundary_poly,
    add = TRUE,
    border = "black",
    lwd = 1,
    col = NA
  )
  
  # ── Full-width colorbar panel ──────────────────────────────────────────
  
  par(
    mar = c(4.0, 6, 0.8, 6),
    bty = "n",
    pty = "m"
  )
  
  pad <- (
    vmax - vmin
  ) * 0.25
  
  cbar_z <- matrix(
    seq_len(n_colors),
    ncol = 1
  )
  
  image(
    x = breaks,
    y = c(0.7, 1.3),
    z = cbar_z,
    col = wind_palette,
    breaks = seq(
      0.5,
      n_colors + 0.5,
      by = 1
    ),
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = c(
      vmin - pad,
      vmax + pad
    ),
    ylim = c(
      0.5,
      1.5
    )
  )
  
  # NOTE: no divider-line loop here -- that's an SLP-specific style
  # (SLP_nicolosi.R uses it, the Python animate_nicolosi_wind.py source
  # this script is translated from does not), left over from the SLP
  # script being used as its starting point. Matches the wind PlateCarree/
  # Robinson scripts too, which also have no colorbar dividers.
  
  rect(
    vmin,
    0.7,
    vmax,
    1.3,
    border = "black",
    lwd = 1.2
  )
  
  axis(
    side = 1,
    at = seq(
      vmin,
      vmax,
      2
    ),
    labels = seq(
      vmin,
      vmax,
      2
    ),
    cex.axis = 1.0
  )
  
  mtext(
    "Wind Speed at 10 m (m/s)",
    side = 1,
    line = 2.7,
    cex = 1.0
  )
  
  # Close the PNG device.
  dev.off()
  
  # Save the completed frame path.
  frame_paths[i] <- fpath
  
  cat(
    sprintf(
      "Saved frame %d/%d\n",
      i,
      n_frames
    )
  )
}

# ── Encode frames as MP4 ─────────────────────────────────────────────────

video_path <- file.path(
  output_dir,
  "nicolosi_wind.mp4"
)

av_encode_video(
  input = frame_paths,
  output = video_path,
  framerate = 4,
  codec = "libx264"
)

cat(
  sprintf(
    "Done! Animation saved to: %s\n",
    video_path
  )
)