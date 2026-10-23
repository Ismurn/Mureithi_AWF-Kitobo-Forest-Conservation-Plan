# 00_setup.R  —  Kitobo Forest Conservation Plan
# Packages, configuration, CRS, directory scaffold, and shared helpers.
# Source this FIRST (the phase scripts do it themselves).
# Paths are relative to the repo root: open the .Rproj in RStudio, or run
# Rscript from the repo root.

# 0a. Packages 
# Work from the repo root (the RStudio project already does this).
if (!file.exists("scripts/00_setup.R") && requireNamespace("here", quietly = TRUE))
  setwd(here::here())

# Non-interactive runs (Rscript) have no CRAN mirror set.
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@"))
  options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  sf, terra, tidyverse, tmap,     # core spatial + data
  exactextractr,                  # zonal statistics
  elevatr,                        # open DEM access
  whitebox,                       # hydrological analysis (file-based)
  trend,                          # Theil-Sen (sens.slope) for recovery trend
  scico, viridis, patchwork,      # visualisation
  units, here                     # areas with units; project-root paths
)

# WhiteboxTools: install the binary ONCE (uncomment on first run), then init each session.
# whitebox::install_whitebox()
wbt_init()

tmap_mode("view")   # interactive maps with satellite basemaps

# 0b. Project CRS 
crs_utm <- "EPSG:32737"   # WGS84 / UTM Zone 37S. Metric WORKING CRS: area, buffers, hydrology, GEE.
crs_pts <- "EPSG:21037"   # CRS the survey CSV is recorded in: Arc 1960 / UTM Zone 37S (Survey of Kenya).
                          # D1 REVISED: read as WGS84 (32737) the points sit ~315 m off the Google Earth
                          # outline. 01_boundary.R re-checks this on every run; see docs/methodology.md s1.
crs_geo <- "EPSG:4326"    # geographic (for FIRMS AOI strings and *_wgs84 layers)

# 0c. EXTERNAL data roots (EDIT THESE to match your machine) 
# Kept as variables so no absolute paths are hardcoded in the pipeline.
dir_tsv   <- "C:/Workspace/Mureithi_RPI_Conservation_Governance/data/Tsavo_Conservation_Areas"  # RanchesTTWCA.shp, SettledAreas.shp
dir_kenya <- "C:/Workspace/Mureithi_RPI_Conservation_Governance/data/Country_Borders"
# dir_kenya holds: kenya_admin_boundary_wgs84.shp, kenya_terrestrial_protected_areas_wgs84.shp

# 0d. Repo-internal paths (relative to repo root) 
path_pts     <- "data/raw/Kitobo_Forest_Coordinates.csv"   # survey points, recorded in crs_pts
path_kml_ref <- "data/raw/Final_Kitobo_Forest.kml"         # Google Earth outline (WGS84): QC reference only, gitignored
dir_proc     <- "data/processed"      # intermediates written by the pipeline
dir_gee      <- "data/raw/gee"        # drop GEE Drive exports (GeoTIFFs) here
dir_shp      <- "outputs/shapefiles"
dir_fig      <- "outputs/figures"
dir_ras      <- "outputs/rasters"
dir_tab      <- "outputs/tables"
dir_wbt      <- "data/processed/wbt"  # WhiteboxTools scratch (inputs/outputs on disk)

for (d in c(dir_proc, dir_gee, dir_shp, dir_fig, dir_ras, dir_tab, dir_wbt))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

#  0e. Tunable analysis parameters
params <- list(
  boundary         = "hull",  # D2 (lead's decision, 2026-10-23): "hull" = concave hull of the survey points;
                              # "outline" = polygon traced in Google Earth (176.5 ha; needs path_kml_ref), kept for comparison. See 01_boundary.R.
  hull_ratio       = 0.10,   # concave-hull ratio for the forest polygon (D2); also drives the sensitivity table.
  catchment_buffer = 5000,   # metres to buffer the forest for the hydrology extent (D4).
  dem_zoom         = 13,     # elevatr tile zoom. z = 13 gives 9.5 m cells here, interpolated from ~30 m
                             # SRTM-class data, so the real resolution is ~30 m.
  stream_area_km2  = 1,      # contributing area (km2) a cell needs to start a stream (D3). TUNE vs imagery.
                             # The old 1000-cell default (0.09 km2) gave 2.9 km/km2: far too dense for this flat, porous plain.
  fire_date        = as.Date("2026-02-15"),  # mid-Feb 2026 (lead's note in the GEE script: pre-fire window ends 10 Feb,
                             # post-fire window starts 1 Mar). Monthly composites only need the right month; refine from FIRMS.
  gee_pre          = c("2026-01-10", "2026-02-10"),   # Sentinel-2 dNBR windows as run in Earth Engine on 2026-10-23
  gee_post         = c("2026-03-01", "2026-04-15")    # (set in the GEE script; recorded here for the results note)
)

# 0f. Secrets 
# FIRMS MAP_KEY is read from .Renviron (gitignored): add a line  FIRMS_MAP_KEY=xxxx
firms_map_key <- Sys.getenv("FIRMS_MAP_KEY")   # empty "" until set; checked in 03_...

# 0g. Helpers 
# Total area in hectares (numeric) of a projected sf/sfc polygon set (0 for an empty set).
area_ha <- function(x) if (length(x) == 0) 0 else sum(as.numeric(units::set_units(sf::st_area(x), "ha")))

# Standard satellite basemap set for tmap view mode (first = default).
basemaps <- c("Esri.WorldImagery", "OpenStreetMap",
              "Esri.WorldGrayCanvas", "Esri.WorldTopoMap")

# Save an interactive tmap as HTML. A self-contained file needs pandoc (bundled with
# RStudio); plain Rscript has none, so fall back to a file plus a _files folder.
save_map <- function(m, file) {
  grDevices::pdf(NULL)            # tmap opens a graphics device while saving: keep it off the disk (no Rplots.pdf)
  on.exit(grDevices::dev.off())
  tmap_save(m, file, selfcontained = rmarkdown::pandoc_available())
}

# Save a base-graphics plot (terra::plot etc.) as a PNG in outputs/figures:
#   save_png("name.png", { plot(...); lines(...) })
# The plotting code runs inside the open device. (ggplot objects use ggsave directly.)
save_png <- function(name, code, width = 8, height = 6.5, res = 300) {
  png(file.path(dir_fig, name), width = width, height = height, units = "in", res = res)
  on.exit(dev.off())
  force(code)
  invisible(file.path(dir_fig, name))
}

message("Setup complete | working CRS = ", crs_utm, " | survey points CRS = ", crs_pts,
        " | FIRMS key ", ifelse(nzchar(firms_map_key), "loaded", "NOT set (needed in 03_...)"))
