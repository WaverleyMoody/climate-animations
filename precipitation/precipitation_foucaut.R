# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: precipitation_foucaut.R
# Description: Generates the precipitation climatology animation from ERA5 reanalysis data (1979-2000), rendered in the Foucaut projection. Weekly totals (cm water equivalent) computed server-side by the derived-era5-single-levels-daily-statistics CDS dataset, averaged into a 22-year climatological mean per 48-frame weekly cycle.
#
# Note: For the Plate Carrée, Robinson, and Nicolosi projections, see the other scripts in the precipitation scripts folder.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Volumes/CLIMATEDATA/precip/climatology/climatology_precip_48frame.nc"
output_dir <- "/Users/waverleymoody/Downloads/precipitation_animation/r_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

proj_name  <- "foucaut"
target_crs <- "+proj=fouc +datum=WGS84 +units=m +no_defs"

frames_dir <- file.path(output_dir, proj_name, "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

# ── Plot settings (identical to the PlateCarree/Robinson precipitation
# scripts) ──────────────────────────────────────────────────────────────
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
# since colorRampPalette() only supports evenly-spaced stops -- needed
# here (unlike the temperature Foucaut script) because precipitation's
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

n_colors       <- 60  # matches the Python script's 60 color steps
precip_palette <- make_custom_palette(precip_colors, precip_positions, n_colors)
breaks         <- seq(vmin, vmax, length.out = n_colors + 1)

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
r_stack <- rast(data_file, subds = "tp_weekly_total_cm")  # already in cm, no conversion needed
r_stack <- rotate(r_stack)              # ERA5 is 0-360 -- fix before reprojecting
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
  
  field   <- r_stack[[i]]
  field_p <- mask(project(field, target_crs), boundary_poly)
  
  label <- frame_labels[i]
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(fpath, width = 1800, height = 980, res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1.3))
  
  # ── Map panel ──
  # bty = "n" reliably suppresses the plot-region box; box=FALSE isn't a
  # real terra::plot() parameter and gets silently dropped.
  par(mar = c(2, 3, 7.5, 2), bty = "n")
  plot(field_p, col = precip_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "")
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.0, col = "black")
  plot(st_geometry(borders_p), add = TRUE, lwd = 0.8, border = "black")
  plot(st_geometry(states_p), add = TRUE, lwd = 0.5, border = "black")
  
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
  
  mtext("ERA-5 | Climate Reanalyzer", side = 3, line = 6, adj = 0,
        cex = 1.0, font = 2)
  mtext(sprintf("%s; 1979\u20132000 Weekly Total Mean", label), side = 3, line = 6,
        adj = 1, cex = 1.0, font = 2)
  
  # ── Colorbar panel ──
  # Custom-drawn via image() rather than terra's plot(..., legend.only=TRUE):
  # that legend mechanism positions itself in whole-device NDC coordinates,
  # not panel-local ones, so it doesn't respect this layout panel. Drawing
  # it directly with image() gives full control and keeps it inside this
  # panel.
  par(mar = c(4.5, 6, 1.2, 6), bty = "n")
  cbar_mat <- matrix(seq(vmin, vmax, length.out = n_colors), ncol = 1)
  # xlim wider than the actual vmin/vmax data range lSeaves blank space on
  # both sides within the panel, so the bar reads as shorter/narrower
  # instead of stretching edge-to-edge across the full frame width.
  #
  # x = breaks (not seq(vmin, vmax, length.out=n_colors)) is important:
  # passing n_colors points makes image() treat them as bin CENTERS, which
  # auto-extends the painted area half a bin-width past both vmin and
  # vmax -- that's what let the purple bleed past the vmax=21 border in
  # the first render. Passing the n_colors+1 breaks instead makes x the
  # true cell edges, so the painted region stops exactly at vmax.
  pad <- (vmax - vmin) * 0.25
  image(x = breaks, y = 1, z = cbar_mat,
        col = precip_palette, breaks = breaks,
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(vmin - pad, vmax + pad))
  usr <- par("usr")
  rect(vmin, usr[3], vmax, usr[4], border = "black", lwd = 1.2)
  axis(1, at = seq(vmin, vmax, tick_step), labels = seq(vmin, vmax, tick_step),
       cex.axis = 1.0)
  mtext("Precipitation (cm)", side = 1, line = 3.0, cex = 1.0)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_path <- file.path(output_dir, "precipitation_climatology_foucaut.mp4")
av_encode_video(frame_paths, output = video_path, framerate = 4, codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))