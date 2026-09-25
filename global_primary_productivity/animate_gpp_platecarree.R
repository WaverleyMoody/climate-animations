# SDSU Climate Informatics Lab
# San Diego State University
# by Waverley Moody
# Supervised by Distinguished Professor Samuel Shen
# R Code Version 1.0.0
#
# A reproduction of the University of Washington General Circulation
# Animations Library by Professor John Michael Wallace.
#
# Script: animate_gpp_platecarree.R
# Description: Generates the Gross Primary Productivity climatology
#     animation from MOD17A2HGF (2000-2009), rendered in the Plate Carree
#     projection.
# Note: For the Robinson, Foucaut, and Nicolosi projections, see the other
#     scripts in the gpp scripts folder. Data excludes 2000-01-01 through
#     2000-02-10 (early Terra commissioning gap) - the animation effectively
#     starts 2000-02-18, 454 frames total. Each frame is written to its own
#     PNG and already-rendered frames are skipped on rerun, so the render can
#     be interrupted and resumed; frames are stitched with av once all are
#     present, then the frame folder is deleted.

library(ncdf4)
library(sf)
library(rnaturalearth)
library(av)

# ---- Config -------------------------------------------------------------
DATA_PATH   <- "/Volumes/CLIMATEDATA/gpp_2000_2009.nc"
OUTPUT_PATH <- "/Volumes/CLIMATEDATA/gpp_2000_2009_platecarree_R.mp4"
FRAMES_DIR  <- "/Volumes/CLIMATEDATA/_frames_tmp_gpp_platecarree_R"
VAR_NAME    <- "gpp"

FPS <- 4.7  # matches the original UW clip: ~454 frames / 98 seconds

# Test mode: render only the first MAX_FRAMES frames into a short "_test" MP4.
# Set to NULL for the full render (test frames are reused).
MAX_FRAMES <- NULL

VMIN  <- 0
VMAX  <- 0.12  # matches the observed value range from spot-checking
GAMMA <- 0.6   # <1 darkens mid-range values, pulling more of the colormap's
# saturated greens into the visible range instead of clustering
# everything in the pale end

# matplotlib "YlGn" (ColorBrewer 9-class, linearly interpolated to 256 colors).
N_COLORS <- 256
PALETTE <- colorRampPalette(c(
  "#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679",
  "#41ab5d", "#238443", "#006837", "#004529"
))(N_COLORS)

# PowerNorm equivalent: colors are evenly spaced in normalized space
# u = ((v - VMIN) / (VMAX - VMIN))^GAMMA, so the data-space breaks are
# VMIN + (VMAX - VMIN) * u^(1 / GAMMA). N_COLORS + 1 true cell edges, not
# midpoints, so no color overflows past VMIN/VMAX (lab colorbar principle).
BREAKS <- VMIN + (VMAX - VMIN) * seq(0, 1, length.out = N_COLORS + 1)^(1 / GAMMA)

OCEAN_COLOR     <- "#0a1a3c"  # dark navy, matches the original clip's ocean
LAND_BASE_COLOR <- "#f2f2ea"  # off-white/gray base for non-vegetated land
COASTLINE_COLOR <- "#808080"  # matplotlib "gray" (R's "gray" is much lighter)

# Figure: 12 x 6 in at 150 dpi, same as the Python figure (even dims for H.264).
FIG_W_PX <- 1800
FIG_H_PX <- 900
FIG_DPI  <- 150
PT_X <- (FIG_DPI / 72) / FIG_W_PX  # one typographic point, in figure fractions
PT_Y <- (FIG_DPI / 72) / FIG_H_PX

# GPP colorbar placement/style. Figure-fraction [left, bottom, width, height].
COLORBAR_RECT       <- c(0.28, 0.085, 0.44, 0.035)
COLORBAR_LABEL      <- expression("GPP (kg C/m"^2 * "/8-day)")
COLORBAR_TICKS      <- c(0, 0.03, 0.06, 0.09, 0.12)
COLORBAR_TEXT_COLOR <- "black"

# Map axes placement, figure-fraction [left, bottom, width, height] -- same
# as the Python MAP_RECT. Plate Carree keeps a 2:1 aspect, so (like Cartopy)
# the map is height-limited and centered horizontally inside this box.
MAP_RECT <- c(0.01, 0.13, 0.98, 0.80)
.map_box <- function() {
  h_in <- MAP_RECT[4] * FIG_H_PX
  w_px <- h_in * 2
  left <- MAP_RECT[1] + (MAP_RECT[3] - w_px / FIG_W_PX) / 2
  c(left, left + w_px / FIG_W_PX, MAP_RECT[2], MAP_RECT[2] + MAP_RECT[4])  # par("plt") order
}

# R's lwd = 1 is 1/96 in (0.75 pt); convert matplotlib point widths.
.pt_lwd <- function(pt) pt / 0.75


# ---- Data loading -------------------------------------------------------

.read_time <- function(nc) {
  vals  <- ncvar_get(nc, "time")
  units <- ncatt_get(nc, "time", "units")$value
  parts <- regmatches(units, regexec("^\\s*(\\w+)\\s+since\\s+(.+?)\\s*$", units))[[1]]
  if (length(parts) != 3) stop("Unrecognized time units: ", units)
  secs <- c(second = 1, seconds = 1, minute = 60, minutes = 60,
            hour = 3600, hours = 3600, day = 86400, days = 86400)[tolower(parts[2])]
  if (is.na(secs)) stop("Unsupported time step in units: ", units)
  origin <- as.POSIXct(sub("T", " ", sub("Z$", "", parts[3])), tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M", "%Y-%m-%d"))
  as.Date(origin + vals * secs)
}

# Reads one time step as an [x, y] matrix -- only a single frame is ever in
# memory (the R equivalent of chunks={"time": 1} + cache_frame_data=False).
.read_slice <- function(nc, t_idx) {
  dn <- vapply(nc$var[[VAR_NAME]]$dim, function(d) d$name, character(1))
  start <- rep(1, length(dn))
  count <- rep(-1, length(dn))
  start[dn == "time"] <- t_idx
  count[dn == "time"] <- 1
  m <- ncvar_get(nc, VAR_NAME, start = start, count = count, collapse_degen = TRUE)
  if (which(dn == "x") > which(dn == "y")) m <- t(m)
  m
}

# Project convention: force a real data read, not just a metadata open.
.verify_readable <- function(nc) {
  if (is.null(nc$var[[VAR_NAME]])) stop("Variable ", VAR_NAME, " missing from ", DATA_PATH)
  m <- .read_slice(nc, 1)
  if (all(is.na(m))) stop("First ", VAR_NAME, " slice in ", DATA_PATH, " is entirely NA")
  rm(m)
  invisible(gc())
}

# Works out how to put the grid in the increasing-x, increasing-y, -180..180
# order that image() needs, once, so each frame just applies it.
.grid_layout <- function(x, y) {
  wrap <- max(x) > 180
  x_std <- if (wrap) ((x + 180) %% 360) - 180 else x
  x_ord <- order(x_std)
  y_ord <- order(y)
  xs <- x_std[x_ord]
  ys <- y[y_ord]
  edges <- function(v) {
    d <- diff(v)
    c(v[1] - d[1] / 2, v[-1] - d / 2, v[length(v)] + d[length(d)] / 2)
  }
  list(x_ord = x_ord, y_ord = y_ord, x_edges = edges(xs), y_edges = edges(ys))
}

.format_date <- function(d) {
  # "Feb 18, 2000" -- month.abb is locale-independent, and as.integer drops the
  # leading zero portably (the Python %-d code is not portable).
  sprintf("%s %d, %s", month.abb[as.integer(format(d, "%m"))],
          as.integer(format(d, "%d")), format(d, "%Y"))
}

.list_frames <- function() {
  # ^ anchor also excludes AppleDouble "._" sidecars on the exFAT drive
  sort(list.files(FRAMES_DIR, pattern = "^frame_\\d{5}\\.png$", full.names = TRUE))
}


# ---- Rendering ----------------------------------------------------------

.draw_gpp_colorbar <- function() {
  # Built directly from PALETTE / VMIN / VMAX / GAMMA, so it can't drift out of
  # sync with the map coloring. Like matplotlib's colorbar for a PowerNorm, the
  # colors are evenly spaced along the bar and the ticks sit at their
  # power-scaled positions.
  l <- COLORBAR_RECT[1]; b <- COLORBAR_RECT[2]
  w <- COLORBAR_RECT[3]; h <- COLORBAR_RECT[4]
  xpos <- function(v) l + ((v - VMIN) / (VMAX - VMIN))^GAMMA * w
  
  edges <- seq(l, l + w, length.out = N_COLORS + 1)
  rect(edges[-(N_COLORS + 1)], b, edges[-1], b + h, col = PALETTE, border = NA)
  rect(l, b, l + w, b + h, border = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.6))
  
  segments(xpos(COLORBAR_TICKS), b, xpos(COLORBAR_TICKS), b - 3.5 * PT_Y,
           col = COLORBAR_TEXT_COLOR, lwd = .pt_lwd(0.8))
  
  tick_label_top <- b - (3.5 + 3.5) * PT_Y   # tick length + matplotlib's default tick pad
  text(xpos(COLORBAR_TICKS), tick_label_top, labels = sprintf("%.2f", COLORBAR_TICKS),
       adj = c(0.5, 1), cex = 9 / 12, col = COLORBAR_TEXT_COLOR)
  text(l + w / 2, tick_label_top - (9 * 1.2 + 6) * PT_Y, labels = COLORBAR_LABEL,
       adj = c(0.5, 1), cex = 11 / 12, col = COLORBAR_TEXT_COLOR)
}

.render_frame <- function(frame_path, z, grid, date_label, land, coast) {
  # Clamp to [VMIN, VMAX] so out-of-range values take the end colors (matching
  # matplotlib's colormap under/over behavior); NA stays transparent so the
  # ocean/land base layers show through.
  z <- z[grid$x_ord, grid$y_ord]
  z <- pmin(pmax(z, VMIN), VMAX)
  
  # Write to a partial file and rename on success, so an interrupted save can
  # never leave a truncated frame that the resume logic would then skip.
  tmp_path <- paste0(frame_path, ".partial")
  png(tmp_path, width = FIG_W_PX, height = FIG_H_PX, res = FIG_DPI, pointsize = 12, bg = "white")
  
  # Map panel, in lon/lat coordinates.
  box <- .map_box()
  par(plt = box, xpd = FALSE)   # plt only -- setting mar here would override it
  plot.new()
  plot.window(xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i")
  
  rect(-180, -90, 180, 90, col = OCEAN_COLOR, border = NA)
  plot(sf::st_geometry(land), col = LAND_BASE_COLOR, border = NA, add = TRUE)
  image(x = grid$x_edges, y = grid$y_edges, z = z, col = PALETTE, breaks = BREAKS,
        add = TRUE, useRaster = TRUE)
  plot(sf::st_geometry(coast), col = COASTLINE_COLOR, lwd = .pt_lwd(0.3), add = TRUE)
  rect(-180, -90, 180, 90, border = "black", lwd = .pt_lwd(0.8))
  
  # Figure-fraction overlay for the titles and colorbar.
  par(plt = c(0, 1, 0, 1), new = TRUE, xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  title_y <- box[4] + 6 * PT_Y   # matplotlib's default title pad (6 pt)
  text(box[1], title_y, "Gross Primary Productivity", adj = c(0, 0), cex = 1)
  text(box[2], title_y, date_label, adj = c(1, 0), cex = 1)
  
  .draw_gpp_colorbar()
  dev.off()
  file.rename(tmp_path, frame_path)
}

.assemble_video <- function(output_path, n_frames = NULL) {
  frame_files <- .list_frames()
  if (!is.null(n_frames)) frame_files <- head(frame_files, n_frames)
  message(sprintf("Assembling %d frames into %s at %s fps...", length(frame_files), output_path, FPS))
  if (file.exists(output_path)) file.remove(output_path)
  av::av_encode_video(frame_files, output = output_path, framerate = FPS,
                      codec = "libx264", verbose = FALSE)
  message("Saved animation to ", output_path)
}


# ---- Main ---------------------------------------------------------------

main <- function() {
  dir.create(FRAMES_DIR, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(FRAMES_DIR, pattern = "\\.partial$", full.names = TRUE))
  
  nc <- nc_open(DATA_PATH)
  on.exit(nc_close(nc), add = TRUE)
  .verify_readable(nc)
  
  dates    <- .read_time(nc)
  n_frames <- length(dates)
  grid     <- .grid_layout(ncvar_get(nc, "x"), ncvar_get(nc, "y"))
  
  # Cartopy's default OCEAN/LAND/COASTLINE features are Natural Earth 110m.
  land  <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
  
  test_mode     <- !is.null(MAX_FRAMES) && MAX_FRAMES < n_frames
  render_frames <- if (test_mode) MAX_FRAMES else n_frames
  if (test_mode) {
    message(sprintf("TEST MODE: rendering the first %d of %d frames.", render_frames, n_frames))
  } else {
    message(sprintf("Rendering %d frames at %s fps...", n_frames, FPS))
  }
  
  already_done <- length(.list_frames())
  if (already_done > 0) {
    message(sprintf("Resuming: %d frame(s) already rendered, skipping those.", already_done))
  }
  
  pb <- txtProgressBar(min = 0, max = render_frames, style = 3)
  for (i in seq_len(render_frames)) {
    setTxtProgressBar(pb, i)
    frame_path <- file.path(FRAMES_DIR, sprintf("frame_%05d.png", i - 1))
    if (file.exists(frame_path)) next
    
    z <- .read_slice(nc, i)
    .render_frame(frame_path, z, grid, .format_date(dates[i]), land, coast)
    rm(z)
    if (i %% 25 == 0) invisible(gc())
  }
  close(pb)
  
  if (test_mode) {
    .assemble_video(sub("\\.mp4$", "_test.mp4", OUTPUT_PATH), n_frames = render_frames)
    message("Set MAX_FRAMES <- NULL and rerun for the full render (test frames are reused).")
    return(invisible(NULL))
  }
  
  rendered_count <- length(.list_frames())
  if (rendered_count == n_frames) {
    .assemble_video(OUTPUT_PATH)
    unlink(FRAMES_DIR, recursive = TRUE)
  } else {
    message(sprintf(paste0("WARNING: expected %d frames but found %d -- skipping assembly. ",
                           "Rerun this script to fill in the missing frame(s) first."),
                    n_frames, rendered_count))
  }
}

main()