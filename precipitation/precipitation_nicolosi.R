# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: precipitation_nicolosi.R
# Description: Generates the precipitation climatology animation from ERA5 reanalysis data (1979-2000), rendered as a double-hemisphere Nicolosi Globular projection. Weekly totals (cm water equivalent) computed server-side by the derived-era5-single-levels-daily-statistics CDS dataset, averaged into a 22-year climatological mean per 48-frame weekly cycle.
#
# Note: For the Plate Carrée, Robinson, and Foucaut projections, see the other scripts in the precipitation scripts folder.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# Disable S2 spherical geometry for hemisphere clipping.
sf::sf_use_s2(FALSE)

# ── Paths ────────────────────────────────────────────────────────────────

data_file <- "/Volumes/CLIMATEDATA/precip/climatology/climatology_precip_48frame.nc"

output_dir <- paste0(
  "/Users/waverleymoody/Downloads/",
  "precipitation_animation/r_output"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

frames_dir <- file.path(
  output_dir,
  "nicolosi",
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

# ── Plot settings ────────────────────────────────────────────────────────

vmin <- 0
vmax <- 21
tick_step <- 3

# Same gradient as the Python Plate Carrée/Robinson/Foucaut/Nicolosi
# precipitation scripts: white at 0 cm blends into blue over just the
# first 0.3 cm, then continues through blue -> green -> yellow -> orange
# -> red -> purple up to 21 cm. Positions are fractions of vmax, mirroring
# the Python color_stops list exactly.
precip_colors <- c(
  "#ffffff", # white           (0 cm)
  "#d0e8f0", # pale blue       (0.3 cm)
  "#7ab8e0", # light blue      (2 cm)
  "#3060c0", # blue            (3 cm)
  "#00994d", # green           (6 cm)
  "#66cc00", # yellow-green    (8 cm)
  "#ffff00", # yellow          (9 cm)
  "#ff9900", # orange          (10 cm)
  "#ff3300", # red-orange      (11 cm)
  "#990000", # dark red        (12 cm)
  "#660066"  # purple          (21 cm)
)
precip_positions <- c(0, 0.3, 2, 3, 6, 8, 9, 10, 11, 12, 21) / vmax

# Mirrors LinearSegmentedColormap.from_list with (position, color) tuples,
# since colorRampPalette() only supports evenly-spaced stops -- needed
# here (unlike the temperature Nicolosi script) because precipitation's
# color stops are deliberately uneven (white compressed into just the
# first 0.3 cm out of a 21 cm range).
make_custom_palette <- function(colors, positions, n) {
  rgb_mat  <- t(col2rgb(colors)) / 255
  eval_pos <- seq(0, 1, length.out = n)
  r <- approx(positions, rgb_mat[, 1], xout = eval_pos)$y
  g <- approx(positions, rgb_mat[, 2], xout = eval_pos)$y
  b <- approx(positions, rgb_mat[, 3], xout = eval_pos)$y
  rgb(r, g, b)
}

n_colors <- 60  # matches the Python script's 60 color steps

precip_palette <- make_custom_palette(
  precip_colors,
  precip_positions,
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

r_stack <- rast(
  data_file,
  subds = "tp_weekly_total_cm"
)

# Convert ERA5 longitudes from 0-360 to -180-180.
r_stack <- rotate(r_stack)

# Already in cm -- no unit conversion needed (unlike temperature's
# Kelvin -> Celsius step here).

n_frames <- nlyr(r_stack)

lons <- xFromCol(
  r_stack,
  1:ncol(r_stack)
)

lats <- yFromRow(
  r_stack,
  1:nrow(r_stack)
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
  
  list(
    col_index = col_index,
    XY = XY,
    proj_str = proj_str,
    template = template,
    boundary_poly = boundary_poly,
    coastline_p = coastline_p,
    borders_p = borders_p
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

# ── Generate PNG frames ──────────────────────────────────────────────────

frame_paths <- character(
  n_frames
)

for (i in seq_len(n_frames)) {
  
  frame_mat <- as.matrix(
    r_stack[[i]],
    wide = TRUE
  )
  
  west_r <- render_hemisphere_raster(
    frame_mat,
    west
  )
  
  east_r <- render_hemisphere_raster(
    frame_mat,
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
      "%s; 1979\u20132000 Weekly Total Mean",
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
    col = precip_palette,
    breaks = breaks,
    range = c(vmin, vmax),
    axes = FALSE,
    legend = FALSE,
    main = "",
    asp = 1
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
    col = precip_palette,
    breaks = breaks,
    range = c(vmin, vmax),
    axes = FALSE,
    legend = FALSE,
    main = "",
    asp = 1
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
  
  # x = breaks (true cell edges, length n_colors+1) rather than n_colors
  # midpoints -- this is the edge-vs-center fix from the Foucaut script:
  # passing midpoints makes image() auto-extend the painted area half a
  # bin-width past both vmin and vmax, letting the purple end bleed past
  # the border rect() drawn at the true vmax.
  image(
    x = breaks,
    y = c(0.7, 1.3),
    z = cbar_z,
    col = precip_palette,
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
      tick_step
    ),
    labels = seq(
      vmin,
      vmax,
      tick_step
    ),
    cex.axis = 1.0
  )
  
  mtext(
    "Precipitation (cm)",
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
  "precipitation_climatology_nicolosi.mp4"
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