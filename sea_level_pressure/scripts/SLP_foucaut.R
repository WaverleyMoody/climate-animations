# SDSU Climate Informatics
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# San Diego State University
#
# R translation of animate_foucaut_SLP.py
# Sea level pressure climatology animation in the Foucaut projection
# (PROJ's `fouc` -- plain Foucaut, NOT `fouc_s`/Foucaut Sinusoidal).
#
# Built directly on the finalized 2m_temp_foucaut.R template (diamond
# boundary masking, custom-drawn colorbar, US-only states, font sizing all
# carry over unchanged) -- swaps in SLP's variable, unit conversion, color
# scale, and the between-band divider lines that animate_SLP.py has but
# animate_2m_temp.py doesn't.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Users/waverleymoody/Downloads/climate_data_by_variable/slp_climatology.nc"
output_dir <- "/Users/waverleymoody/Downloads/slp_animation_R"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

proj_name  <- "foucaut"
target_crs <- "+proj=fouc +datum=WGS84 +units=m +no_defs"

frames_dir <- file.path(output_dir, proj_name, "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

# ── Plot settings (identical to animate_SLP.py's fixed 970-1040 mb scale) ──
vmin <- 970
vmax <- 1040
color_stops <- c(
  '#2d004b', # dark purple  (970 mb)
  '#800080', # purple
  '#8b008b', # magenta
  '#ff00ff', # bright pink
  '#003153', # dark blue
  '#b0c4de', # blue
  '#ffffff', # white
  '#ff6600', # orange
  '#ff0000', # red
  '#8b0000'  # dark red     (1040 mb)
)
n_colors <- 60
slp_palette <- colorRampPalette(color_stops)(n_colors)
breaks <- seq(vmin, vmax, length.out = n_colors + 1)

# ── Natural Earth basemap layers ───────────────────────────────────────────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
# US only -- see 2m_temp_foucaut.R for why: ne_states() with no country
# filter pulls every country's first-level admin boundaries worldwide.
states    <- ne_states(country = "United States of America", returnclass = "sf")

# Reproject once, not per frame
coastline_p <- st_transform(coastline, target_crs)
borders_p   <- st_transform(borders, target_crs)
states_p    <- st_transform(states, target_crs)

# ── Load pre-computed climatology ──────────────────────────────────────────
r_stack <- rast(data_file, subds = "msl")
r_stack <- rotate(r_stack)              # ERA5 is 0-360 -- fix before reprojecting
r_stack <- r_stack / 100                # Pa -> hPa/mb
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
# See 2m_temp_foucaut.R for the full explanation -- terra::project() fills
# the rectangular bounding box and extrapolates outside the true domain
# instead of returning NA, so we trace the boundary and mask() to it.
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
  par(mar = c(2, 3, 7.5, 2), bty = "n")
  plot(field_p, col = slp_palette, breaks = breaks, range = c(vmin, vmax),
       axes = FALSE, legend = FALSE, main = "")
  
  plot(st_geometry(coastline_p), add = TRUE, lwd = 1.0, col = "black")
  plot(st_geometry(borders_p), add = TRUE, lwd = 0.8, border = "black")
  plot(st_geometry(states_p), add = TRUE, lwd = 0.5, border = "black")
  
  for (m in meridians) plot(m, add = TRUE, col = "gray", lwd = 0.3)
  for (p in parallels) plot(p, add = TRUE, col = "gray", lwd = 0.3)
  
  # Pole labels -- single points, not a longitude loop (same as temp version)
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
  # Same custom image()-drawn approach as the temp version (terra's
  # legend.only positions itself in whole-device coordinates, not
  # panel-local ones, so it doesn't respect this layout panel).
  par(mar = c(4.5, 6, 1.2, 6), bty = "n")
  pad <- (vmax - vmin) * 0.25
  # x is given as the actual cell EDGES (breaks, n_colors+1 values spanning
  # exactly vmin to vmax), not cell centers -- passing centers left it to
  # image() to compute edges as midpoints and auto-extend the two outer
  # edges half a cell-width beyond vmin/vmax, which is what let the color
  # poke past the border rect drawn at the exact vmin/vmax values.
  cbar_z <- matrix(seq_len(n_colors), ncol = 1)
  image(x = breaks, y = c(0.5, 1.5), z = cbar_z,
        col = slp_palette, breaks = seq(0.5, n_colors + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(vmin - pad, vmax + pad))
  
  # Divider lines between the colormap's 60 discrete bands -- present in
  # animate_SLP.py/animate_foucaut_SLP.py but not the temp version, so
  # added here to match. [-c(1, length)] drops the two outer edges, same
  # as Python's np.linspace(0, 1, 61)[1:-1] slice.
  for (b in breaks[-c(1, length(breaks))]) {
    abline(v = b, col = "black", lwd = 0.5)
  }
  
  usr <- par("usr")
  rect(vmin, usr[3], vmax, usr[4], border = "black", lwd = 1.2)
  axis(1, at = seq(vmin, vmax, 10), labels = seq(vmin, vmax, 10), cex.axis = 1.0)
  mtext("Pressure at Mean Sea Level (mb)", side = 1, line = 3.0, cex = 1.0)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_path <- file.path(output_dir, "slp_climatology_foucaut.mp4")
av_encode_video(frame_paths, output = video_path, framerate = 4, codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))