# 01_boundary.R  —  Build the Kitobo Forest polygon
# Reads the survey points (Arc 1960 / UTM 37S), checks the CRS against the
# operator's Google Earth outline, reprojects to WGS84 / UTM 37S, builds the
# boundary polygon (concave hull of the points, or the operator's traced
# outline), puts the forest in its landscape context, and saves the
# authoritative shapefile used by every downstream step.
# Prereq: none (sources 00_setup.R itself).
source("scripts/00_setup.R")

# 1a. Read the survey points 
# Projected X/Y recorded in crs_pts (NOT EPSG:4326). Exact repeats (duplicate GPS
# fixes) are dropped.
Kitobo.Raw <- read_csv(path_pts, show_col_types = FALSE) %>% distinct(X, Y)
message("Boundary points: ", nrow(Kitobo.Raw), " unique")

# 1b. CRS check against the reference outline 
# The Google Earth outline (WGS84) is an independent trace of the same forest. Read in
# the right CRS the points sit on that line (median offset of a few metres); read in
# the wrong one they are hundreds of metres away.
Kitobo.Ref <- NULL
if (file.exists(path_kml_ref)) {
  Kitobo.Ref <- st_read(path_kml_ref, quiet = TRUE) %>% st_zm() %>% st_transform(crs_utm)
  crs_check <- purrr::map_dfr(unique(c(crs_pts, "EPSG:21037", "EPSG:32737")), function(crs_try) {
    d <- Kitobo.Raw %>%
      st_as_sf(coords = c("X", "Y"), crs = crs_try) %>%
      st_transform(crs_utm) %>%
      st_distance(Kitobo.Ref) %>%
      as.numeric()
    tibble(points_read_as = crs_try, median_m = median(d), p90_m = unname(quantile(d, 0.9)))
  })
  print(crs_check)
  readr::write_csv(crs_check, file.path(dir_tab, "01_crs_check.csv"))
  if (crs_check$points_read_as[which.min(crs_check$median_m)] != crs_pts)
    warning("The points line up better with the reference outline in another CRS than ",
            "crs_pts (", crs_pts, "). Check crs_pts in 00_setup.R.", call. = FALSE)
} else {
  message("NOTE: reference outline ", path_kml_ref, " not found; skipping the CRS check.")
}

# 1c. Points in the working CRS 
Kitobo.Pts <- Kitobo.Raw %>%
  st_as_sf(coords = c("X", "Y"), crs = crs_pts) %>%
  st_transform(crs_utm)

# 1d. Boundary candidates
# (i)  D2: CONCAVE hull of the survey points. The points lie only on the perimeter, so a low
#      ratio cuts chords across the interior. Needs sf >= 1.0-14 with GEOS >= 3.11.
# (ii) The operator's Google Earth trace is an ORDERED outline of the same forest: close it,
#      node it where it crosses itself, polygonize, and drop slivers (< 0.1 ha) left by the
#      thin loop at the north tip. No hull needed.
hull_at <- function(r) st_make_valid(st_concave_hull(st_union(Kitobo.Pts), ratio = r))

outline_polygon <- function(ref) {
  xy     <- st_coordinates(ref)[, 1:2]
  closed <- st_sfc(st_linestring(rbind(xy, xy[1, ])), crs = st_crs(ref))
  pieces <- st_cast(st_collection_extract(st_polygonize(st_node(closed)), "POLYGON"), "POLYGON")
  st_make_valid(st_union(pieces[vapply(seq_along(pieces), function(i) area_ha(pieces[i]), 0) > 0.1]))
}
iou <- function(a, b) area_ha(st_intersection(a, b)) / area_ha(st_union(c(a, b)))

Kitobo.Hull    <- hull_at(params$hull_ratio)
Kitobo.Outline <- if (!is.null(Kitobo.Ref)) outline_polygon(Kitobo.Ref)
if (!is.null(Kitobo.Outline))   # kept so 03 can report burned area under both boundaries
  st_write(st_sf(name = "Kitobo traced outline", geometry = Kitobo.Outline),
           file.path(dir_proc, "kitobo_outline.gpkg"), append = FALSE, quiet = TRUE)

boundary <- switch(params$boundary,
  hull    = Kitobo.Hull,
  outline = if (is.null(Kitobo.Outline)) stop("params$boundary = 'outline' needs ", path_kml_ref) else Kitobo.Outline,
  stop("params$boundary must be 'hull' or 'outline'"))
Kitobo.Forest <- st_sf(name = "Kitobo Forest", geometry = boundary)

message("Kitobo Forest area = ", round(area_ha(Kitobo.Forest), 1), " ha (boundary = ", params$boundary,
        if (params$boundary == "hull") paste0(", hull_ratio = ", params$hull_ratio), ")")
if (!is.null(Kitobo.Outline))
  message("The traced outline encloses ", round(area_ha(Kitobo.Outline), 1), " ha; the hull at ratio ",
          params$hull_ratio, " matches it with IoU ", round(iou(Kitobo.Hull, Kitobo.Outline), 2))

# How the hull responds to the ratio and how well it matches the traced outline. The forest
# is usually quoted at ~150 ha; compare that with the traced outline before settling.
ratios   <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00)
hulls    <- purrr::map(ratios, hull_at)
no_ref   <- is.null(Kitobo.Outline)
hull_tab <- tibble(
  hull_ratio             = ratios,
  hull_ha                = purrr::map_dbl(hulls, area_ha),
  iou_vs_outline         = if (no_ref) NA_real_ else purrr::map_dbl(hulls, ~ iou(.x, Kitobo.Outline)),
  outline_missed_ha      = if (no_ref) NA_real_ else purrr::map_dbl(hulls, ~ area_ha(st_difference(Kitobo.Outline, .x))),
  hull_beyond_outline_ha = if (no_ref) NA_real_ else purrr::map_dbl(hulls, ~ area_ha(st_difference(.x, Kitobo.Outline)))
)
print(hull_tab, width = Inf)
readr::write_csv(hull_tab, file.path(dir_tab, "01_hull_ratio_sensitivity.csv"))

# Hulls at four ratios drawn over the area enclosed by the traced outline.
if (!no_ref) {
  sel     <- ratios %in% c(0.10, 0.30, 0.50, 1.00)
  hull_sf <- st_sf(hull = sprintf("ratio %.2f: %.0f ha, IoU %.2f",
                                  ratios[sel], hull_tab$hull_ha[sel], hull_tab$iou_vs_outline[sel]),
                   geometry = do.call(c, hulls[sel]))
  cmp <- ggplot() +
    geom_sf(data = Kitobo.Outline, fill = "goldenrod", alpha = 0.35, colour = NA) +
    geom_sf(data = hull_sf, fill = NA, colour = "darkgreen", linewidth = 0.7) +
    geom_sf(data = Kitobo.Pts, colour = "red", size = 0.5) +
    facet_wrap(~hull) +
    labs(title = "Concave hull of the survey points vs the traced outline",
         subtitle = paste0("Yellow = area enclosed by the Google Earth trace (",
                           round(area_ha(Kitobo.Outline)), " ha); red = survey points")) +
    theme_minimal(base_size = 9) + theme(axis.text = element_blank(), panel.grid = element_blank())
  ggsave(file.path(dir_fig, "01_hull_vs_outline.png"), cmp, width = 9, height = 6, dpi = 200)
}

# 1e. Static check figure (no basemap needed) 
# Survey points under both CRS readings vs the Google Earth outline.
Kitobo.Pts.asWGS84 <- Kitobo.Raw %>% st_as_sf(coords = c("X", "Y"), crs = crs_utm)
lab_ok  <- "Survey points read as Arc 1960 / UTM 37S"
lab_old <- "Survey points read as WGS84 / UTM 37S (old D1)"
lab_ref <- "Google Earth outline (WGS84)"
chk <- ggplot()
if (!no_ref) chk <- chk + geom_sf(data = Kitobo.Outline, fill = "goldenrod", alpha = 0.15, colour = NA)
chk <- chk +
  geom_sf(data = Kitobo.Forest, fill = "forestgreen", alpha = 0.25, colour = "darkgreen") +
  geom_sf(data = Kitobo.Pts.asWGS84, aes(colour = lab_old), shape = 4, size = 1.8) +
  geom_sf(data = Kitobo.Pts, aes(colour = lab_ok), size = 1.8)
if (!no_ref) chk <- chk + geom_sf(data = Kitobo.Ref, aes(colour = lab_ref), linewidth = 0.7)
chk <- chk +
  scale_colour_manual(NULL, values = setNames(c("red", "grey40", "goldenrod"), c(lab_ok, lab_old, lab_ref))) +
  labs(title = "Kitobo Forest: survey points vs Google Earth outline",
       subtitle = paste0("Green = boundary used (", params$boundary,
                         if (params$boundary == "hull") paste0(", ratio ", params$hull_ratio),
                         "); pale yellow = area enclosed by the traced outline")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(dir_fig, "01_boundary_crs_check.png"), chk, width = 9, height = 5.5, dpi = 200)

#  1f. Interactive QC: boundary + points + outline on satellite imagery 
# Inspect this map. If the polygon clips real forest or bridges a gap, change
# params$boundary / params$hull_ratio in 00_setup.R and re-run from the top.
qc_map <- tm_shape(Kitobo.Forest) +
  tm_polygons(fill = "forestgreen", fill_alpha = 0.25, col = "darkgreen", lwd = 2, id = "name")
if (!no_ref) qc_map <- qc_map + tm_shape(Kitobo.Ref) + tm_lines(col = "gold", lwd = 2)
qc_map <- qc_map +
  tm_shape(Kitobo.Pts) + tm_dots(fill = "red", size = 0.4) +
  tm_title("Kitobo Forest: survey points (red), boundary (green), Google Earth outline (gold)") +
  tm_basemap(server = basemaps)
print(qc_map)

# 1g. Landscape context: host ranch / settlement / protected area
# Kitobo lies WEST of the TTWCA ranch and settlement layers, so for those the useful answer
# is the distance to the nearest feature (0 = the forest intersects it).
context_tab <- tryCatch(
  purrr::imap_dfr(
    list("TTWCA ranches"   = list(file.path(dir_tsv,   "RanchesTTWCA.shp"), "R_NAME"),
         "Settled areas"   = list(file.path(dir_tsv,   "SettledAreas.shp"), "R_NAME"),
         "Protected areas" = list(file.path(dir_kenya, "kenya_terrestrial_protected_areas_wgs84.shp"), "NAME")),
    function(src, layer) {
      lay <- st_read(src[[1]], quiet = TRUE) %>% st_make_valid()
      d   <- as.numeric(st_distance(st_transform(Kitobo.Forest, st_crs(lay)), lay))
      tibble(layer = layer, nearest = as.character(lay[[src[[2]]]][which.min(d)]),
             distance_km = round(min(d) / 1000, 1))
    }),
  error = function(e) {
    message("NOTE: landscape layers not loaded (set dir_tsv / dir_kenya in 00_setup.R). Skipping 1g. [",
            conditionMessage(e), "]")
    NULL
  })
if (!is.null(context_tab)) {
  print(context_tab)
  readr::write_csv(context_tab, file.path(dir_tab, "01_landscape_context.csv"))
}

# 1h. Save the authoritative shapefile + a GeoPackage copy
st_write(Kitobo.Forest, file.path(dir_shp, "Kitobo_Forest.shp"),
         append = FALSE, quiet = TRUE)
st_write(Kitobo.Forest, file.path(dir_proc, "kitobo_forest.gpkg"),
         append = FALSE, quiet = TRUE)
save_map(qc_map, file.path(dir_fig, "01_boundary_qc.html"))

message("Saved: ", file.path(dir_shp, "Kitobo_Forest.shp"),
        "  -> upload this to GEE asset projects/ee-imnduta1/assets/Kitobo_Forest")
