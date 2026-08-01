# SDSU Climate Informatics
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# San Diego State University

# A reproduction of the University of Washington General Circulation
# Animations Library, originally created by Professor John Michael Wallace.

# Script: precipitable_water_robinson.R
# Description: Generates the precipitable water (total column water vapour)
# climatology animation from ERA5 reanalysis data (1979-2000), rendered
# in the Robinson projection.

library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(av)

# ── Paths ────────────────────────────────────────────────────────────────
data_file  <- "/Users/waverleymoody/Downloads/climate_data_by_variable/precip_water_climatology.nc"
output_dir <- "/Users/waverleymoody/Downloads/precipitable_water_animation_R"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Projection to render ────────────────────────────────────────────────────
proj_name  <- "robinson"
target_crs <- "+proj=robin +datum=WGS84 +units=m +no_defs"

# ── Plot settings ────────────────────────────────────────────────────────
vmin <- 0
vmax <- 80

# Same 15 color stops as the Python custom_pw colormap, interpolated to 40
color_stops <- c(
  '#d3d3d3', # light gray      (0 kg/m2)
  '#808080', # gray            (6 kg/m2)
  '#4b2922', # brown           (12 kg/m2)
  '#d2b48c', # tan/light brown (20 kg/m2)
  '#adcfe6', # light blue      (24 kg/m2)
  '#30306a', # blue            (30 kg/m2)
  '#004A00', # dark green      (35 kg/m2)
  '#00cc00', # green           (40 kg/m2)
  '#ffff00', # yellow          (50 kg/m2)
  '#ff6600', # orange          (56 kg/m2)
  '#ff0000', # red             (60 kg/m2)
  '#690909', # dark red
  '#67013e', # dark pink
  '#ff00ff', # magenta         (80 kg/m2)
  '#38014b'  # dark purple
)
n_colors <- 40
pw_palette <- colorRampPalette(color_stops)(n_colors)
breaks <- seq(vmin, vmax, length.out = n_colors + 1)

# ── Natural Earth basemap layers (cached locally after first download) ─────
cat("Downloading/loading basemap layers...\n")
coastline <- ne_coastline(scale = 50, returnclass = "sf")
borders   <- ne_countries(scale = 50, returnclass = "sf")
states    <- ne_states(country = "United States of America", returnclass = "sf")

# ── Load pre-computed climatology ──────────────────────────────────────────
r_stack <- rast(data_file, subds = "tcwv")
r_stack <- rotate(r_stack)   # ERA5 uses 0-360 longitude; convert to -180/180
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
  field   <- r_stack[[i]]
  field_p <- project(field, target_crs, method = "bilinear")
  
  label <- as.character(frame_labels[i])
  fpath <- file.path(frames_dir, sprintf("frame_%03d.png", i - 1))
  
  png(fpath, width = 1800, height = 980, units = "px", res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(5, 1))
  
  # ── Panel 1: map ──────────────────────────────────────────────────────
  par(mar = c(2, 3, 7.5, 2))
  
  plot(field_p, col = pw_palette, breaks = breaks, range = c(vmin, vmax),
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
  mtext("ERA-5 | Climate Reanalyzer", side = 3, line = 6, adj = 0,
        cex = 0.9, font = 2)
  mtext(sprintf("%s; 1979\u20132000 Weekly Mean", label), side = 3, line = 6,
        adj = 1, cex = 0.9, font = 2)
  
  # ── Panel 2: custom horizontal colorbar ──────────────────────────────
  par(mar = c(3, 6, 0.5, 6))
  n_swatches <- length(pw_palette)
  plot(NA, xlim = c(0, n_swatches), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rect(0:(n_swatches - 1), 0, 1:n_swatches, 1, col = pw_palette, border = NA)
  box(col = "black", lwd = 1)
  
  tick_vals <- seq(vmin, vmax, 10)
  tick_pos  <- (tick_vals - vmin) / (vmax - vmin) * n_swatches
  axis(1, at = tick_pos, labels = tick_vals, cex.axis = 0.8, tck = -0.3,
       mgp = c(3, 0.5, 0))
  mtext("Precipitable Water (kg m-2)", side = 1, line = 2, cex = 0.9)
  
  dev.off()
  
  frame_paths[i] <- fpath
  cat(sprintf("Saved frame %d/%d\n", i, n_frames))
}

video_name <- "precipitable_water_climatology_robinson.mp4"
video_path <- file.path(output_dir, video_name)

av_encode_video(frame_paths, output = video_path, framerate = 4,
                codec = "libx264")

cat(sprintf("Done! Animation saved to: %s\n", video_path))