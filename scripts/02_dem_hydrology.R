# 02_dem_hydrology.R  —  Terrain + surface hydrology for Kitobo
# DEM over the buffered catchment -> terrain derivatives (terra) -> WhiteboxTools
# flow routing, stream network, channel heads ("river starting points"), and TWI.
# Prereq: 01_boundary.R has produced the forest polygon.
# NOTE: modelled streams are SURFACE drainage, not the springs feeding this
#       groundwater forest — report them as such; foreground TWI for water story.
source("scripts/00_setup.R")

Kitobo.Forest <- st_read(file.path(dir_proc, "kitobo_forest.gpkg"), quiet = TRUE) %>%
  st_transform(crs_utm)

# 2a. Buffered catchment extent (flow originates upslope, D4) 
Kitobo.Buff <- st_buffer(Kitobo.Forest, dist = params$catchment_buffer)

# 2b. DEM acquisition (open terrain tiles)
# elevatr pulls blended open terrain tiles (SRTM-based). For publication, swap in
# FABDEM (canopy/building-bias removed) — see docs/methodology.md — and read it with rast().
Kitobo.DEM <- get_elev_raster(Kitobo.Buff, z = params$dem_zoom, clip = "locations") %>%
  rast() %>%
  project(crs_utm)

dem_path      <- file.path(dir_wbt, "dem.tif")
breached_path <- file.path(dir_wbt, "dem_breached.tif")
writeRaster(Kitobo.DEM, dem_path, overwrite = TRUE)

# 2c. Terrain derivatives (fast, in-memory)
Kitobo.Terrain <- terrain(Kitobo.DEM,
                          v = c("slope", "aspect", "TPI", "TRI"),
                          unit = "degrees")
writeRaster(Kitobo.Terrain, file.path(dir_ras, "kitobo_terrain.tif"), overwrite = TRUE)

# 2d. Hydro-conditioning: breach depressions (preferred over filling) 
wbt_breach_depressions_least_cost(dem = dem_path, output = breached_path, dist = 20)

# 2e. Flow direction (D8 pointer) + flow accumulation (cells) 
ptr_path  <- file.path(dir_wbt, "d8_pointer.tif")
facc_path <- file.path(dir_wbt, "facc_cells.tif")
wbt_d8_pointer(dem = breached_path, output = ptr_path)
wbt_d8_flow_accumulation(input = breached_path, output = facc_path, out_type = "cells")

# 2f. Extract streams -> Strahler order -> vectorise 
# params$stream_area_km2 sets drainage density: the contributing area a cell needs to start a
# channel, converted to cells with the real DEM resolution. TUNE against imagery/known channels.
stream_cells <- round(params$stream_area_km2 * 1e6 / prod(res(Kitobo.DEM)))
message("Stream threshold: ", params$stream_area_km2, " km2 = ", stream_cells, " cells of ",
        round(res(Kitobo.DEM)[1], 1), " m")
streams_ras <- file.path(dir_wbt, "streams.tif")
order_ras   <- file.path(dir_wbt, "strahler.tif")
streams_vec <- file.path(dir_wbt, "streams.shp")
wbt_extract_streams(flow_accum = facc_path, output = streams_ras,
                    threshold = stream_cells)
wbt_strahler_stream_order(d8_pntr = ptr_path, streams = streams_ras, output = order_ras)
wbt_raster_streams_to_vector(streams = streams_ras, d8_pntr = ptr_path, output = streams_vec)

Kitobo.Streams <- st_read(streams_vec, quiet = TRUE) %>% st_set_crs(crs_utm)

# 2g. "River starting points" = channel heads (upstream-most nodes) 
# Pragmatic rule: the first vertex of a segment that is not shared with any other
# vertex in the network is a source node (a channel head). Refine if the network
# has many junctions (e.g. snap-and-count within a small tolerance).
seg_ids    <- seq_len(nrow(Kitobo.Streams))
seg_starts <- Kitobo.Streams %>%
  mutate(.seg = seg_ids) %>%
  { suppressWarnings(st_cast(., "POINT")) } %>%   # "repeating attributes" is expected here
  group_by(.seg) %>% slice(1) %>% ungroup()
all_verts  <- st_cast(st_geometry(Kitobo.Streams), "POINT")
Kitobo.Sources <- seg_starts %>%
  filter(lengths(st_equals(st_geometry(.), all_verts)) == 1)

# Clip the network + heads to the forest for the within-forest view
Kitobo.Streams.F <- suppressWarnings(st_intersection(Kitobo.Streams, Kitobo.Forest))
st_write(Kitobo.Streams, file.path(dir_shp, "Kitobo_Streams.shp"), append = FALSE, quiet = TRUE)
st_write(Kitobo.Sources, file.path(dir_shp, "Kitobo_ChannelHeads.shp"), append = FALSE, quiet = TRUE)

# 2h. TWI (saturation index) — needs specific catchment area + slope
sca_path   <- file.path(dir_wbt, "sca.tif")
slope_path <- file.path(dir_wbt, "slope.tif")
twi_path   <- file.path(dir_ras, "kitobo_twi.tif")
wbt_d8_flow_accumulation(input = breached_path, output = sca_path,
                         out_type = "specific contributing area")
wbt_slope(dem = breached_path, output = slope_path, units = "degrees")
wbt_wetness_index(sca = sca_path, slope = slope_path, output = twi_path)
Kitobo.TWI <- rast(twi_path)

# 2i. Map: TWI + modelled streams + channel heads over the forest 
hydro_map <- tm_shape(crop(Kitobo.TWI, vect(Kitobo.Buff))) +
  tm_raster(col.scale = tm_scale_intervals(values = "viridis", style = "quantile", n = 7),
            col_alpha = 0.85,
            col.legend = tm_legend(title = "TWI (wetness)")) +
  tm_shape(Kitobo.Streams) + tm_lines(col = "steelblue", lwd = 1.3) +
  tm_shape(Kitobo.Sources) + tm_dots(fill = "red", size = 0.3) +
  tm_shape(Kitobo.Forest)  + tm_borders(col = "darkgreen", lwd = 2.5) +
  tm_title("Kitobo — modelled surface drainage, channel heads & wetness") +
  tm_basemap(server = basemaps)
print(hydro_map)
save_map(hydro_map, file.path(dir_fig, "02_hydrology.html"))

# 2j. Static map for reports (no basemap): hillshade, streams, heads, forest
png(file.path(dir_fig, "02_hydrology.png"), width = 1500, height = 1250, res = 150)
hs <- shade(terrain(Kitobo.DEM, "slope", unit = "radians"),
            terrain(Kitobo.DEM, "aspect", unit = "radians"), angle = 40, direction = 270)
plot(hs, col = gray.colors(100, 0.25, 0.95), legend = FALSE,
     main = "Kitobo: modelled surface drainage (blue), channel heads (red), forest (green)")
plot(Kitobo.DEM, add = TRUE, alpha = 0.35, col = terrain.colors(60))
lines(vect(Kitobo.Streams), col = "steelblue", lwd = 0.7)
points(vect(Kitobo.Sources), col = "red", pch = 16, cex = 0.4)
lines(vect(Kitobo.Forest), col = "darkgreen", lwd = 2.5)
invisible(dev.off())

# 2k. Drainage summary: use it to tune params$stream_area_km2
dem_km2   <- as.numeric(global(!is.na(Kitobo.DEM), "sum")) * prod(res(Kitobo.DEM)) / 1e6
stream_km <- sum(as.numeric(st_length(Kitobo.Streams))) / 1000
message("Streams: ", nrow(Kitobo.Streams), " segments, ", round(stream_km, 1), " km, drainage density ",
        round(stream_km / dem_km2, 2), " km/km2 over ", round(dem_km2, 1), " km2")
message("Hydrology done. Channel heads: ", nrow(Kitobo.Sources), " (",
        sum(lengths(st_intersects(Kitobo.Sources, Kitobo.Forest)) > 0), " inside the forest)",
        " | tune params$stream_area_km2 if drainage density looks wrong.")
