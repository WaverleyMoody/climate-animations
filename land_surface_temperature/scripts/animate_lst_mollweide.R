# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: animate_lst_mollweide.R
#
# Description: Generates the Land Surface Temperature (Day) animation from NASA Terra/MODIS
# MOD11C3 v6.1 (Feb 2000-Dec 2013, native 0.05 deg CMG grid), rendered in the
# Mollweide projection. Unlike the other ERA5-based variables in this library,
# this is a 167-frame real monthly time series, not a 48-frame climatological composite.
#
# Note: Mollweide has a working forward AND inverse transform via PROJ, so unlike Nicolosi
# this does not require the manual rasterize/focal gap-fill pipeline -- terra::project()
# on a per-frame raster is sufficient (mirrors the Cartopy GeoAxes path used in the Python
# version, which is why this script follows the same standard rendering approach as Foucaut).

library(terra)
library(ncdf4)
library(sf)
library(rnaturalearth)
library(av)

IN_PATH  <- "/Volumes/CLIMATEDATA/modis_lst_monthly/lst_day_timeseries_2000_2013_native05deg_celsius.nc"
OUT_DIR  <- "/Volumes/CLIMATEDATA/modis_lst_monthly/animation_output/"
OUT_MP4  <- file.path(OUT_DIR, "lst_day_mollweide_2000_2013_R.mp4")
FRAME_DIR <- file.path(OUT_DIR, "_frames_tmp_R")

VAR_NAME <- "LST_Day"
VMIN <- -50
VMAX <- 55   # degrees C, physically sane global land LST range (converted from 220-330 K)

# Custom colormap: light blue -> dark blue -> purple -> bright purple -> red -> orange -> yellow
CUSTOM_COLORS <- c(
  "#ADD8E6",  # light blue
  "#00008B",  # dark blue
  "#800080",  # purple
  "#DA00FF",  # bright purple
  "#FF0000",  # red
  "#FFA500",  # orange
  "#FFFF00"   # yellow
)
N_COLORS <- 256
PALETTE  <- colorRampPalette(CUSTOM_COLORS)(N_COLORS)
# Cell edges (length N_COLORS + 1), not midpoints -- required so image() does not let
# colors overflow past the VMIN/VMAX border (see lab colorbar principle).
BREAKS <- seq(VMIN, VMAX, length.out = N_COLORS + 1)

CEAN_GRAY <- "#D3D3D3"  # ocean / no-data background color

FPS <- 6                 # slower than a 48-frame climatology loop since this has 167 real frames
DPI <- 150
FIG_W_IN <- 12
FIG_H_IN <- 7

MOLL_CRS <- "+proj=moll +datum=WGS84 +units=m +no_defs"


# Confirms the input time-series file is a real, uncorrupted read before animating --
# project convention. Forces an actual ncvar_get() read, not just a shallow file open.
.verify_readable <- function(path, varname) {
  nc <- tryCatch(nc_open(path), error = function(e) e)
  if (inherits(nc, "error")) {
    return(list(ok = FALSE, reason = conditionMessage(nc)))
  }
  on.exit(nc_close(nc), add = TRUE)
  result <- tryCatch({
    v <- ncvar_get(nc, varname, start = c(1, 1, 1), count = c(1, 1, 1))
    list(ok = TRUE, reason = NA_character_)
  }, error = function(e) list(ok = FALSE, reason = conditionMessage(e)))
  result
}


.format_date <- function(time_posix) {
  # Time values are stored as month starts (day = 1); format to "April 2000" style.
  format(time_posix, "%B %Y")
}


# Builds a Mollweide-projected mask raster: value 1 inside the valid global map ellipse,
# NA everywhere outside it. Used to fill the map interior gray (matching the Python
# ax.set_facecolor(gray) behavior, which only shows within the clipped GeoAxes boundary)
# without needing manual ellipse polygon geometry.
.build_ocean_mask <- function(template_r, moll_crs) {
  ones_r <- template_r
  values(ones_r) <- 1
  crs(ones_r) <- "EPSG:4326"
  project(ones_r, moll_crs, method = "near")
}


.render_frame <- function(i, r_stack, times, moll_crs, mask_moll,
                          coast_moll, borders_moll, frame_dir) {
  frame_native <- r_stack[[i]]
  crs(frame_native) <- "EPSG:4326"
  frame_moll <- project(frame_native, moll_crs, method = "bilinear")
  
  png_path <- file.path(frame_dir, sprintf("frame_%04d.png", i))
  png(png_path, width = FIG_W_IN * DPI, height = FIG_H_IN * DPI, res = DPI, bg = "white")
  
  # Layout: main map on top, thin horizontal colorbar strip below (mirrors
  # orientation="horizontal", shrink=0.7 in the Python cbar).
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.88, 0.12))
  
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot(ext(mask_moll), col = NA, border = NA, axes = FALSE, xlab = "", ylab = "")
  
  # Gray fill for the map interior (ocean + any masked/no-data land pixels), then the
  # data raster on top with NA cells transparent, so gray shows through gaps.
  plot(mask_moll, col = OCEAN_GRAY, colNA = NA, legend = FALSE, axes = FALSE, add = TRUE)
  plot(frame_moll, col = PALETTE, breaks = BREAKS, colNA = NA,
       legend = FALSE, axes = FALSE, add = TRUE)
  
  plot(st_geometry(coast_moll), add = TRUE, col = "black", lwd = 0.5)
  plot(st_geometry(borders_moll), add = TRUE,
       col = adjustcolor("black", alpha.f = 0.5), lwd = 0.3)
  
  mtext("MODIS Land Surface Temperature\n(Daytime)", side = 3, line = -2, adj = 0,
        font = 2, cex = 0.85)
  mtext(.format_date(times[i]), side = 3, line = -2, adj = 1, font = 2, cex = 0.85)
  
  # Horizontal colorbar strip built from BREAKS as cell edges, per lab convention.
par(mar = c(3.5, 4, 0.5, 4))
image(x = BREAKS, y = 1, z = matrix(BREAKS[-length(BREAKS)], ncol = 1),
      col = PALETTE, breaks = BREAKS, axes = FALSE, xlab = "", ylab = "")
axis(1, at = pretty(BREAKS), cex.axis = 0.8, line = 0)
mtext(expression("Land Surface Temperature (" * degree * "C)"), side = 1, line = 2.3, cex = 0.75)
  
  dev.off()
  invisible(png_path)
}


main <- function() {
  if (!file.exists(IN_PATH)) {
    stop(sprintf("ERROR: input time-series file not found: %s\nRun concat_lst_timeseries.R first.", IN_PATH))
  }
  
  verify <- .verify_readable(IN_PATH, VAR_NAME)
  if (!verify$ok) {
    stop(sprintf("ERROR: input file failed verification read: %s", verify$reason))
  }
  
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FRAME_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # terra::rast on a NetCDF is file-backed / lazy by default -- frames are pulled into
  # memory one at a time inside .render_frame(), mirroring the chunks={"time":1} approach
  # used in the Python/xarray pipeline for memory safety on the 8 GB MacBook Air.
  r_stack <- rast(IN_PATH, subds = VAR_NAME)
  n_frames <- nlyr(r_stack)
  times <- time(r_stack)
  cat(sprintf("Loaded %d frames, %s to %s\n", n_frames,
              format(min(times), "%Y-%m-%d"), format(max(times), "%Y-%m-%d")))
  
  ocean_mask <- .build_ocean_mask(r_stack[[1]], MOLL_CRS)
  
  coast_sf <- ne_coastline(scale = 110, returnclass = "sf")
  borders_sf <- ne_countries(scale = 110, returnclass = "sf")
  coast_moll <- st_transform(coast_sf, MOLL_CRS)
  borders_moll <- st_transform(borders_sf, MOLL_CRS)
  
  cat(sprintf("Rendering %d frames...\n", n_frames))
  frame_paths <- character(n_frames)
  for (i in seq_len(n_frames)) {
    frame_paths[i] <- .render_frame(
      i, r_stack, times, MOLL_CRS, ocean_mask, coast_moll, borders_moll, FRAME_DIR
    )
    if (i %% 20 == 0 || i == n_frames) {
      cat(sprintf("  frame %d / %d\n", i, n_frames))
    }
  }
  
  av_encode_video(frame_paths, output = OUT_MP4, framerate = FPS,
                  vfilter = "scale=trunc(iw/2)*2:trunc(ih/2)*2")
  
  unlink(FRAME_DIR, recursive = TRUE)
  
  cat(sprintf("\nSaved animation to %s\n", OUT_MP4))
  
  if (file.exists(OUT_MP4) && file.info(OUT_MP4)$size > 0) {
    cat("PASS -- output file exists and is non-empty.\n")
  } else {
    stop("FAIL -- output file missing or empty.")
  }
}


main()