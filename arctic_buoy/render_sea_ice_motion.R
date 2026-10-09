# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: render_sea_ice_motion.R
#
# Description: Renders the "Sea ice motion" animation - daily Arctic sea ice concentration
#              shading (NOAA/NSIDC G02202 V6) with IABP buoy tracks (red circles, 30-day
#              green fading tails) drawn on top, in the sea ice grid's native NSIDC Polar
#              Stereographic North projection. R translation of render_sea_ice_motion.py.
#
# Data: All inputs are provided in the GitHub Release `cryosphere_v1.0_data`:
#         - sea_ice_concentration_1979_1989.zip, _1990_1999.zip, _2000_2009.zip, _2010_2022.zip
#         - sea_ice_concentration_filled.zip
#         - iabp_daily_positions_1979_2022_projected.parquet
#       Unzip everything into one folder and set DATA_ROOT below to that folder. Expected layout:
#         DATA_ROOT/
#             sea_ice_concentration/<year>/sic_psn25_YYYYMMDD_*.nc
#             sea_ice_concentration_filled/<year>/sic_psn25_YYYYMMDD_*.nc
#             iabp_buoys/iabp_daily_positions_1979_2022_projected.parquet
#       The finished video is written to DATA_ROOT/sea_ice_motion/.
#
# Packages: ncdf4, sf, rnaturalearth, arrow (parquet reader), av
#           install.packages(c("ncdf4", "sf", "rnaturalearth", "arrow", "av"))
#
# Note: TEST_MODE renders a single year (2020) to validate colors, tail rendering, and
#       frame timing before committing to the full 1979-2022 range (~16,000 frames).
#       Set TEST_MODE <- FALSE and adjust YEAR_RANGE to run the full animation.
#
# Note: The whole scene is rotated 90 degrees counterclockwise from the sea ice grid's
#       native orientation. The ice matrix is rotated by index reordering (exact, no
#       resampling), buoy positions by (x, y) -> (-y, x), and all lon/lat layers
#       (coastlines, land, gridlines, labels) are projected through a copy of the grid CRS
#       with lon_0 shifted by -90 degrees - for a polar stereographic projection that shift
#       is an exact 90-degree CCW rotation of the plane, so every layer stays registered.
#
# Note: Days where the source record has no valid ice concentration (SMMR/SSM/I outages,
#       e.g. Jul-Aug 1984) are replaced by the linearly interpolated files in
#       sea_ice_concentration_filled/. Those frames are labeled on screen.
#
# Note: Frames are written as PNGs to local disk and encoded with av at the end (the av
#       package has no streaming writer). The full run needs roughly 5 GB of free local
#       disk for frames; they are deleted after the video is verified and moved. With
#       RESUME <- TRUE, a crashed run picks up from the last finished frame.
#       Run with: caffeinate -i Rscript render_sea_ice_motion.R

suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
  library(rnaturalearth)
  library(arrow)
  library(av)
})
sf_use_s2(FALSE)  # planar lon/lat cropping before projecting; avoids s2 edge cases

TEST_MODE  <- FALSE
YEAR_RANGE <- if (TEST_MODE) 2020:2020 else 1979:2022
TAIL_DAYS  <- 30
RESUME     <- TRUE  # skip frames already rendered by a previous (crashed) run

# --- Paths: set DATA_ROOT to the folder where you unzipped the release assets ---
DATA_ROOT <- "/Volumes/CLIMATEDATA"

ICE_DIR      <- file.path(DATA_ROOT, "sea_ice_concentration")
FILLED_DIR   <- file.path(DATA_ROOT, "sea_ice_concentration_filled")
BUOY_PARQUET <- file.path(DATA_ROOT, "iabp_buoys", "iabp_daily_positions_1979_2022_projected.parquet")

LOCAL_STAGING <- file.path(path.expand("~"), "CLIMATE ANIMATIONS", "staging", "sea_ice_motion")  # local disk scratch
FRAME_DIR     <- file.path(LOCAL_STAGING, if (TEST_MODE) "_frames_tmp_R_test" else "_frames_tmp_R")
FINAL_DEST    <- file.path(DATA_ROOT, "sea_ice_motion")
OUT_NAME      <- if (TEST_MODE) "sea_ice_motion_test_2020_R.mp4" else "sea_ice_motion_1979_2022_R.mp4"

FPS      <- 52   # ~365 daily frames / 7 seconds per year
FIG_W_IN <- 10   # landscape - the grid's bounding box is wider than tall once rotated
FIG_H_IN <- 7
DPI      <- 120
PT       <- function(pt) pt / 12  # matplotlib point sizes -> cex (png pointsize = 12)

TAIL_COLOR_RGB <- c(57, 255, 20)  # "#39FF14" bright/neon green, visible against blue-white ice
TAIL_LWD       <- 1.6             # ~ matplotlib linewidth 1.2 pt
BUOY_CEX       <- 0.45            # ~ matplotlib markersize 3
INTERP_LABEL   <- "Ice interpolated (satellite data gap)"
LAND_GRAY      <- "#dcdcdc"       # land / missing shows as neutral gray, not white (=100% ice)
GRID_GRAY      <- "#808080"       # matplotlib "gray" (R's "gray" is much lighter)

# matplotlib Blues_r (ColorBrewer Blues, reversed), 256 levels over 0-1 concentration.
BLUES_R  <- colorRampPalette(rev(c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6",
                                   "#4292c6", "#2171b5", "#08519c", "#08306b")))(256)
# Cell edges (length n + 1), not midpoints, so colors never overflow past vmin/vmax.
BREAKS   <- seq(0, 1, length.out = length(BLUES_R) + 1)

# Representative points for on-map geographic labels (lon, lat).
LABELS <- data.frame(
  name = c("ALASKA", "GREENLAND", "NORWAY", "SIBERIA"),
  lon  = c(-155, -42, 11, 105),
  lat  = c(67, 72, 64, 72)
)
PARALLEL_LABEL_LON <- -135  # meridian the inline 60/70/80N labels sit on (straight down after rotation)


# Build a CRS identical to `crs` except with its central meridian (lon_0) shifted by -90
# degrees. For polar stereographic this is an exact 90-degree CCW rotation of the plane
# about the pole, so lon/lat layers projected through it land exactly at (x, y) -> (-y, x)
# of their native-grid positions, matching the rotated ice and buoys.
rotate_crs_90ccw <- function(crs) {
  p4  <- crs$proj4string
  pat <- "\\+lon_0=-?[0-9]+\\.?[0-9]*"
  m   <- regmatches(p4, regexpr(pat, p4))
  if (length(m) == 0) stop("Could not find +lon_0= in this CRS's proj4 string: ", p4)
  old_lon0 <- as.numeric(sub("\\+lon_0=", "", m))
  st_crs(sub(pat, paste0("+lon_0=", old_lon0 - 90), p4))
}

# Prefer the gap-filled copy of a day's ice file if one exists in FILLED_DIR.
resolve_ice_path <- function(path) {
  filled <- file.path(FILLED_DIR, basename(dirname(path)), basename(path))
  if (file.exists(filled) && !startsWith(basename(filled), "._")) {
    list(path = filled, filled = TRUE)
  } else {
    list(path = path, filled = FALSE)
  }
}

get_ice_var_name <- function(nc) {
  vars <- names(nc$var)
  for (cand in c("cdr_seaice_conc", "seaice_conc_cdr", "nsidc_nt_seaice_conc")) {
    if (cand %in% vars) return(cand)
  }
  # fall back to the first 2D+ data variable that isn't a flag/qc field
  for (v in vars) {
    if (!grepl("qc|flag", tolower(v)) && nc$var[[v]]$ndims >= 2) return(v)
  }
  stop("Could not identify the ice concentration variable among: ", paste(vars, collapse = ", "))
}

# Returns list(values = [x, y] matrix, x, y). ncvar_get applies scale_factor and masks
# _FillValue, matching xarray's default decoding in the Python version.
load_ice_frame <- function(path, cache) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  if (is.null(cache$var_name)) cache$var_name <- get_ice_var_name(nc)
  vn <- cache$var_name
  
  v <- ncvar_get(nc, vn, collapse_degen = TRUE)  # drops the length-1 time dim
  dn <- vapply(nc$var[[vn]]$dim, function(d) d$name, "")
  dn <- dn[dn %in% c("x", "y")]
  if (identical(dn[1], "y")) v <- t(v)  # guarantee [x, y] ordering
  
  # Defensive rescale: CDR concentration should be a 0-1 fraction. If the file instead
  # stores 0-100 percent (varies by product version), rescale so BREAKS stay correct.
  finite <- v[is.finite(v)]
  if (length(finite) && max(finite) > 1.5) v <- v / 100
  
  list(values = v, x = as.numeric(nc$dim$x$vals), y = as.numeric(nc$dim$y$vals))
}

read_grid_crs <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  for (cand in c("crs", "polar_stereographic", "projection", "Polar_Stereographic_Grid")) {
    for (att in c("crs_wkt", "spatial_ref")) {
      a <- tryCatch(ncatt_get(nc, cand, att), error = function(e) list(hasatt = FALSE))
      if (isTRUE(a$hasatt)) return(st_crs(a$value))
    }
  }
  stop("No grid-mapping variable with crs_wkt/spatial_ref found in ", path)
}

# Rotate the [x, y] value matrix 90 deg CCW into image orientation (rows top->bottom,
# cols left->right) and map to colors. Rotated: image rows follow x descending, image
# columns follow y descending - pure index reordering, the data itself is untouched.
ice_to_raster <- function(values, x, y) {
  v <- values[order(x, decreasing = TRUE), order(y, decreasing = TRUE), drop = FALSE]
  v <- pmin(pmax(v, 0), 1)  # matplotlib clips under/over to the end colors
  cols <- BLUES_R[findInterval(v, BREAKS, all.inside = TRUE)]
  cols[is.na(v)] <- LAND_GRAY
  as.raster(matrix(cols, nrow = nrow(v)))
}

fmt_lon <- function(lo) {
  lo <- ((lo + 180) %% 360) - 180
  if (lo == 0) return("0\u00b0")
  if (abs(lo) == 180) return("180\u00b0")
  if (lo > 0) paste0(lo, "\u00b0E") else paste0(-lo, "\u00b0W")
}

# Text with a white stroke, like matplotlib's patheffects.withStroke.
halo_text <- function(x, y, labels, cex, r, font = 2, col = "black", bg = "white") {
  for (t in seq(0, 2 * pi, length.out = 17)[-17]) {
    text(x + r * cos(t), y + r * sin(t), labels, cex = cex, font = font, col = bg)
  }
  text(x, y, labels, cex = cex, font = font, col = col)
}

# Confirms the encoded video is real and complete before it leaves local disk --
# project convention (verify at every write checkpoint).
.verify_video <- function(path, expected_frames) {
  if (!file.exists(path) || file.size(path) == 0) stop("Video missing or empty: ", path)
  info <- av_media_info(path)
  n <- info$video$frames
  if (is.null(n) || is.na(n) || n < expected_frames) {
    stop(sprintf("Video has %s frames, expected %d: %s", n, expected_frames, path))
  }
  message(sprintf("  verified %s: %d frames, %dx%d, %.1f s",
                  basename(path), n, info$video$width, info$video$height, info$duration))
}


main <- function() {
  if (!dir.exists(ICE_DIR)) {
    stop(ICE_DIR, " not found - set DATA_ROOT at the top of this script to the folder ",
         "where you unzipped the cryosphere_v1.0_data release assets.")
  }
  if (!file.exists(BUOY_PARQUET)) {
    stop(BUOY_PARQUET, " not found - place iabp_daily_positions_1979_2022_projected.parquet in ",
         dirname(BUOY_PARQUET), ".")
  }
  dir.create(FRAME_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FINAL_DEST, recursive = TRUE, showWarnings = FALSE)
  
  message("Collecting daily sea ice files for the requested year range...")
  ice_files <- character(0)
  for (year in YEAR_RANGE) {
    year_dir <- file.path(ICE_DIR, year)
    if (!dir.exists(year_dir)) {
      message("  [WARN] ", year_dir, " does not exist, skipping ", year)
      next
    }
    f <- sort(list.files(year_dir, pattern = "^sic_psn25_.*_v.*\\.nc$", full.names = TRUE))
    ice_files <- c(ice_files, f[!startsWith(basename(f), "._")])
  }
  if (length(ice_files) == 0) stop("No sea ice files found for ", paste(range(YEAR_RANGE), collapse = "-"), " under ", ICE_DIR)
  resolved <- lapply(ice_files, resolve_ice_path)
  message(sprintf("  %d daily ice files found (%d replaced by gap-filled copies)",
                  length(ice_files), sum(vapply(resolved, `[[`, TRUE, "filled"))))
  if (!dir.exists(FILLED_DIR)) {
    message("  [WARN] ", FILLED_DIR, " not found - unzip sea_ice_concentration_filled.zip ",
            "into DATA_ROOT, or gap days render blank")
  }
  
  message("Loading buoy positions from ", BUOY_PARQUET, " ...")
  buoys  <- as.data.frame(read_parquet(BUOY_PARQUET, col_select = c("BuoyID", "date", "x_m", "y_m")))
  b_day  <- as.integer(as.Date(buoys$date))
  o      <- order(b_day, buoys$BuoyID)  # date-sorted so each frame's window is a contiguous slice
  b_day  <- b_day[o]
  b_id   <- buoys$BuoyID[o]
  # Rotate buoy positions 90 degrees CCW to match the rotated ice grid:
  # (x, y) -> (-y, x) is the exact rotation matrix for +90 degrees counterclockwise.
  b_xd   <- -buoys$y_m[o]
  b_yd   <-  buoys$x_m[o]
  message(sprintf("  %d buoy-day rows loaded", nrow(buoys)))
  rm(buoys)
  
  # --- Grid geometry, CRS, and the rotated frame (computed once) ---
  cache <- new.env()
  g0 <- load_ice_frame(ice_files[1], cache)
  x <- g0$x; y <- g0$y
  h <- abs(x[2] - x[1]) / 2  # half cell, for raster edges
  grid_crs <- read_grid_crs(ice_files[1])
  proj <- rotate_crs_90ccw(grid_crs)
  message("  using projection: ", grid_crs$Name, " (rotated 90 deg CCW)")
  
  # Same framing as Python's rotated_extent (cell centers): (-ymax, -ymin, xmin, xmax).
  xlim <- c(-max(y), -min(y))
  ylim <- c(min(x), max(x))
  dx <- diff(xlim); dy <- diff(ylim)
  ras_edges <- c(-(max(y) + h), min(x) - h, -(min(y) - h), max(x) + h)  # xl, yb, xr, yt
  
  # Fit an equal-aspect map box inside the region Python reserves with
  # subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.03), centered like Cartopy.
  avail <- c(l = 0.03, r = 0.97, b = 0.03, t = 0.88)
  aw <- (avail[["r"]] - avail[["l"]]) * FIG_W_IN
  ah <- (avail[["t"]] - avail[["b"]]) * FIG_H_IN
  if (aw / ah > dx / dy) { bh <- ah; bw <- ah * dx / dy } else { bw <- aw; bh <- aw * dy / dx }
  cx <- mean(avail[c("l", "r")]) * FIG_W_IN
  cy <- mean(avail[c("b", "t")]) * FIG_H_IN
  PLT <- c((cx - bw / 2) / FIG_W_IN, (cx + bw / 2) / FIG_W_IN,
           (cy - bh / 2) / FIG_H_IN, (cy + bh / 2) / FIG_H_IN)
  units_per_in <- dx / bw
  halo_r <- (1.25 / 72) * units_per_in  # 2.5 pt stroke width -> 1.25 pt radius
  
  # --- Static vector layers, projected once (Cartopy's default 110m Natural Earth) ---
  crop_nh <- function(g) st_crop(st_make_valid(g), xmin = -180, xmax = 180, ymin = 20, ymax = 90)
  land  <- suppressWarnings(st_transform(crop_nh(st_geometry(ne_countries(scale = 110, returnclass = "sf"))), proj))
  coast <- suppressWarnings(st_transform(crop_nh(st_geometry(ne_coastline(scale = 110, returnclass = "sf"))), proj))
  
  mer_lons <- seq(-180, 150, by = 30)
  meridians <- st_transform(st_sfc(lapply(mer_lons, function(lo) st_linestring(cbind(lo, seq(20, 90, by = 0.5)))), crs = 4326), proj)
  parallels <- st_transform(st_sfc(lapply(c(60, 70, 80), function(la) st_linestring(cbind(seq(-180, 180, by = 0.5), la))), crs = 4326), proj)
  gridlines <- c(meridians, parallels)
  
  # Meridian labels where each meridian crosses the map border (Cartopy x_inline=False).
  box_line <- st_cast(st_as_sfc(st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]), crs = proj)), "LINESTRING")
  pad <- 0.012 * dy
  mer_lab <- do.call(rbind, lapply(seq_along(mer_lons), function(i) {
    p <- suppressWarnings(st_intersection(meridians[i], box_line))
    if (length(p) == 0 || st_is_empty(p)) return(NULL)
    xy <- st_coordinates(p)[, 1:2, drop = FALSE]
    do.call(rbind, lapply(seq_len(nrow(xy)), function(k) {
      px <- xy[k, 1]; py <- xy[k, 2]
      side <- which.min(abs(c(px - xlim[1], xlim[2] - px, py - ylim[1], ylim[2] - py)))
      off  <- list(c(-pad, 0, 1, 0.5), c(pad, 0, 0, 0.5), c(0, -pad, 0.5, 1), c(0, pad, 0.5, 0))[[side]]
      data.frame(x = px + off[1], y = py + off[2], ax = off[3], ay = off[4], lab = fmt_lon(mer_lons[i]))
    }))
  }))
  
  # Parallel labels inline along one meridian (Cartopy y_inline=True).
  par_pts <- st_coordinates(st_transform(st_sfc(lapply(c(60, 70, 80), function(la) st_point(c(PARALLEL_LABEL_LON, la))), crs = 4326), proj))
  par_lab <- data.frame(x = par_pts[, 1], y = par_pts[, 2], lab = paste0(c(60, 70, 80), "\u00b0N"))
  par_lab <- par_lab[par_lab$x > xlim[1] & par_lab$x < xlim[2] & par_lab$y > ylim[1] & par_lab$y < ylim[2], ]
  
  geo_lab <- st_coordinates(st_transform(st_as_sf(LABELS, coords = c("lon", "lat"), crs = 4326), proj))
  
  # --- Render one frame to PNG ---
  render_frame <- function(png_path, ice, frame_date, is_filled) {
    png(png_path, width = FIG_W_IN * DPI, height = FIG_H_IN * DPI, res = DPI, pointsize = 12, bg = "white")
    on.exit(dev.off())
    par(plt = PLT, xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim, ylim)
    
    rasterImage(ice_to_raster(ice$values, ice$x, ice$y),
                ras_edges[1], ras_edges[2], ras_edges[3], ras_edges[4], interpolate = FALSE)
    plot(land, col = LAND_GRAY, border = NA, add = TRUE)
    plot(coast, col = "black", lwd = 0.7, add = TRUE)
    plot(gridlines, col = adjustcolor(GRID_GRAY, 0.6), lty = 2, lwd = 0.8, add = TRUE)
    text(par_lab$x, par_lab$y, par_lab$lab, cex = PT(7), col = GRID_GRAY)
    
    # Buoys in the 30-day window: contiguous slice of the date-sorted arrays.
    fd <- as.integer(frame_date)
    lo <- findInterval(fd - TAIL_DAYS, b_day) + 1  # first row with date > window_start
    hi <- findInterval(fd, b_day)                  # last row with date <= frame_date
    if (hi > lo) {
      k  <- lo:hi
      k  <- k[order(b_id[k], b_day[k])]
      id <- b_id[k]; xs <- b_xd[k]; ys <- b_yd[k]
      r  <- rle(id)
      n_row <- rep(r$lengths, r$lengths)
      pos   <- sequence(r$lengths) - 1L
      m  <- length(id)
      s  <- which(id[-m] == id[-1])  # segment starts (consecutive rows of the same buoy)
      if (length(s)) {
        # Fading tail: oldest segment most transparent, most recent most opaque.
        a <- 0.15 + 0.7 * (pos[s] / pmax(n_row[s] - 2, 1))
        segments(xs[s], ys[s], xs[s + 1], ys[s + 1], lwd = TAIL_LWD,
                 col = rgb(TAIL_COLOR_RGB[1], TAIL_COLOR_RGB[2], TAIL_COLOR_RGB[3],
                           alpha = round(a * 255), maxColorValue = 255))
      }
      last <- cumsum(r$lengths)[r$lengths >= 2]  # current position, buoys with >= 2 days only
      points(xs[last], ys[last], pch = 21, bg = "red", col = "darkred", cex = BUOY_CEX, lwd = 0.6)
    }
    
    halo_text(geo_lab[, 1], geo_lab[, 2], LABELS$name, cex = PT(9), r = halo_r)
    box(lwd = 0.8)
    
    # Outside the map (xpd = NA), at the same axes-fraction heights as the Python version.
    for (j in seq_len(nrow(mer_lab))) {
      text(mer_lab$x[j], mer_lab$y[j], mer_lab$lab[j], adj = c(mer_lab$ax[j], mer_lab$ay[j]),
           cex = PT(7), col = GRID_GRAY, xpd = NA)
    }
    text(xlim[1], ylim[2] + 0.05 * dy, "International Arctic Buoy Program",
         adj = c(0, 0), cex = PT(13), font = 2, xpd = NA)
    text(xlim[2], ylim[2] + 0.05 * dy,
         sprintf("%d/%d/%d", as.integer(format(frame_date, "%m")),
                 as.integer(format(frame_date, "%d")), as.integer(format(frame_date, "%Y"))),
         adj = c(1, 0), cex = PT(13), xpd = NA)
    if (is_filled) {
      # Small gray note under the date, still above the map, so viewers know the ice
      # field on this frame is interpolated rather than observed.
      text(xlim[2], ylim[2] + 0.015 * dy, INTERP_LABEL, adj = c(1, 0), cex = PT(8),
           col = "gray40", xpd = NA)
    }
  }
  
  # --- Frame loop ---
  frames <- file.path(FRAME_DIR, sprintf("frame_%05d.png", seq_along(ice_files)))
  n_skipped <- 0L
  for (i in seq_along(ice_files)) {
    if (RESUME && file.exists(frames[i]) && file.size(frames[i]) > 0) {
      n_skipped <- n_skipped + 1L
      next
    }
    frame_date <- as.Date(strsplit(basename(ice_files[i]), "_")[[1]][3], "%Y%m%d")  # sic_psn25_{YYYYMMDD}_...
    ice <- load_ice_frame(resolved[[i]]$path, cache)
    render_frame(frames[i], ice, frame_date, resolved[[i]]$filled)
    
    if (i %% 30 == 0 || i == length(ice_files)) {
      message(sprintf("  rendered %d/%d frames", i, length(ice_files)))
    }
  }
  if (n_skipped) message(sprintf("  resumed: %d frames already on disk were reused", n_skipped))
  
  # --- Encode on local disk, verify, then move to the exFAT drive ---
  message("Encoding video...")
  local_out <- file.path(LOCAL_STAGING, OUT_NAME)
  av_encode_video(frames, output = local_out, framerate = FPS, verbose = FALSE)
  .verify_video(local_out, length(frames))
  
  # file.rename() fails across filesystems (APFS -> exFAT); copy, confirm, then delete.
  final_out <- file.path(FINAL_DEST, OUT_NAME)
  if (!file.copy(local_out, final_out, overwrite = TRUE) ||
      file.size(final_out) != file.size(local_out)) {
    stop("Copy to ", final_out, " failed or is incomplete; local copy kept at ", local_out)
  }
  unlink(local_out)
  unlink(FRAME_DIR, recursive = TRUE)  # only after the video is verified and safely moved
  message("Done: ", final_out)
}

main()