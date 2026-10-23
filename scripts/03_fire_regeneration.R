# =============================================================================
# 03_fire_regeneration.R  —  Fire extent/severity + regeneration monitoring
# Part A: FIRMS sweep to find the 2026 fire date(s).
# Part B: import Sentinel-2 dNBR/severity (from scripts/gee/01_kitobo_dNBR.js) and map.
# Part C: regeneration — NBR/NDVI recovery series + Theil-Sen trend.
# Prereq: 01_boundary.R (polygon). Parts B/C need the GEE exports in data/raw/gee/.
# =============================================================================
source("scripts/00_setup.R")

Kitobo.Forest <- st_read(file.path(dir_proc, "kitobo_forest.gpkg"), quiet = TRUE) %>%
  st_transform(crs_utm)
Kitobo.Buff <- st_buffer(Kitobo.Forest, dist = params$catchment_buffer)

# =============================================================================
# PART A — FIRMS: find WHEN it burned (375 m VIIRS: good for date, not extent)
# =============================================================================
if (!nzchar(firms_map_key)) {
  message("FIRMS_MAP_KEY not set — skipping Part A. Add it to .Renviron and re-run.")
} else {
  # AOI string "west,south,east,north" in lon/lat, from the forest's own bbox
  bb  <- st_bbox(st_transform(Kitobo.Buff, crs_geo))
  aoi <- sprintf("%.4f,%.4f,%.4f,%.4f", bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])

  # FIRMS area API: DAY_RANGE is 1..5 (checked against https://firms.modaps.eosdis.nasa.gov/api/area/,
  # NOT 10), so sweep the year in 5-day steps. _SP = standard processing (archive, lags the present
  # by months); _NRT = near real time (recent months). Ask for both and de-duplicate.
  starts  <- seq(as.Date("2026-01-01"), Sys.Date(), by = 5)
  sources <- c("VIIRS_SNPP_SP", "VIIRS_NOAA20_SP", "VIIRS_NOAA21_SP",
               "VIIRS_SNPP_NRT", "VIIRS_NOAA20_NRT", "VIIRS_NOAA21_NRT")

  # One call -> list(data, msg). Errors come back as HTTP 400 with a plain-text body (a bad key
  # gives "Invalid MAP_KEY."), so fetch the status and body ourselves instead of letting
  # read_csv() swallow them and report "no fires".
  firms_get <- function(src, d) {
    url <- sprintf("https://firms.modaps.eosdis.nasa.gov/api/area/csv/%s/%s/%s/5/%s",
                   firms_map_key, src, aoi, format(d))
    r <- tryCatch(curl::curl_fetch_memory(url), error = function(e) e)
    if (inherits(r, "error"))
      return(list(data = NULL, msg = paste("request failed:", conditionMessage(r))))
    txt <- strsplit(rawToChar(r$content), "\r?\n")[[1]]
    if (r$status_code != 200 || (length(txt) && !grepl("^latitude,", txt[1])))
      return(list(data = NULL, msg = paste0("HTTP ", r$status_code, ": ", paste(utils::head(txt, 1), collapse = ""))))
    if (length(txt) < 2) return(list(data = NULL, msg = NA_character_))   # empty or header only: no detections
    list(data = dplyr::mutate(readr::read_csv(I(txt), show_col_types = FALSE), source = src),
         msg  = NA_character_)
  }
  calls <- purrr::pmap(tidyr::expand_grid(src = sources, d = starts), firms_get)
  fires <- dplyr::bind_rows(purrr::map(calls, "data"))
  if (nrow(fires))   # SP and NRT overlap for the same satellite: keep one row per detection
    fires <- dplyr::distinct(fires, dplyr::across(dplyr::any_of(
      c("latitude", "longitude", "acq_date", "acq_time", "satellite"))), .keep_all = TRUE)
  msgs <- table(stats::na.omit(purrr::map_chr(calls, "msg")))
  if (length(msgs)) { message("FIRMS answered with text instead of data (count, message):"); print(msgs) }

  if (nrow(fires)) {
    print(dplyr::arrange(dplyr::count(fires, acq_date), dplyr::desc(n)))  # read fire window here
    Kitobo.Fires <- st_as_sf(fires, coords = c("longitude", "latitude"), crs = crs_geo) %>%
      st_transform(crs_utm) %>% st_filter(Kitobo.Buff)
    st_write(Kitobo.Fires, file.path(dir_proc, "kitobo_firms.gpkg"),
             append = FALSE, quiet = TRUE)
    message("FIRMS detections near Kitobo: ", nrow(Kitobo.Fires),
            " — set the GEE pre/post windows from the dates above.")
  } else {
    message("No FIRMS detections returned for 2026 in this AOI (try NRT sources / wider AOI).")
  }
}

# =============================================================================
# PART B — Sentinel-2 dNBR: WHERE & HOW SEVERE (import GEE export, then map)
# Run scripts/gee/01_kitobo_dNBR.js first; download the GeoTIFFs into data/raw/gee/.
# =============================================================================
dnbr_tif <- file.path(dir_gee, "Kitobo_dNBR.tif")
sev_tif  <- file.path(dir_gee, "Kitobo_Severity.tif")

if (file.exists(sev_tif)) {
  Kitobo.dNBR <- rast(dnbr_tif)
  Kitobo.Sev  <- rast(sev_tif)

  sev_map <- tm_shape(Kitobo.Sev) +
    tm_raster(col.scale = tm_scale_categorical(
        values = c("0" = "#d9d9d9", "1" = "#ffffb2", "2" = "#fecc5c",
                   "3" = "#fd8d3c", "4" = "#e31a1c"),
        labels = c("Unburned", "Low", "Mod-low", "Mod-high", "High")),
      col.legend = tm_legend(title = "Burn severity (dNBR)")) +
    tm_shape(Kitobo.Forest) + tm_borders(col = "black", lwd = 2.5) +
    tm_title("Kitobo Forest — 2026 burn severity (Sentinel-2 dNBR)") +
    tm_basemap(server = basemaps)
  print(sev_map)
  save_map(sev_map, file.path(dir_fig, "03_burn_severity.html"))

  # Burned-area fraction by severity class, within the forest (exactextractr idiom)
  burn_frac <- exact_extract(Kitobo.Sev,
                             st_transform(Kitobo.Forest, crs(Kitobo.Sev)),
                             fun = "frac")
  readr::write_csv(as.data.frame(burn_frac), file.path(dir_tab, "burn_fraction_by_class.csv"))
  print(burn_frac)
} else {
  message("Part B skipped: ", sev_tif, " not found. Run scripts/gee/01 and download exports to data/raw/gee/.")
}

# =============================================================================
# PART C — Regeneration: NBR/NDVI recovery over the burn area + Theil-Sen trend
# Run scripts/gee/02_kitobo_regeneration.js to export a monthly NBR stack to data/raw/gee/.
# =============================================================================
nbr_stack_tif <- file.path(dir_gee, "Kitobo_NBR_monthly.tif")   # bands = months, named YYYY_MM_dd_NBR

if (file.exists(nbr_stack_tif)) {
  NBR.Stack <- rast(nbr_stack_tif)
  NBR.Stack[NBR.Stack == 0] <- NA   # GEE can write masked (cloudy) pixels as 0; real NBR is never exactly 0

  # Mean NBR over the forest per month -> a recovery time series
  nbr_series <- exact_extract(NBR.Stack,
                              st_transform(Kitobo.Forest, crs(NBR.Stack)),
                              fun = "mean") %>%
    as.numeric()
  dates <- as.Date(substr(names(NBR.Stack), 1, 10), format = "%Y_%m_%d")   # band names start with the month
  if (all(is.na(dates)))
    stop("Band names are not dates: re-export with scripts/gee/02 (names must start YYYY_MM_dd).")
  recovery <- tibble(date = dates, nbr = nbr_series) %>% arrange(date) %>% drop_na()

  # Theil-Sen slope on POST-fire NBR = a robust regeneration-rate estimate.
  # Set fire_date from Part A before interpreting.
  # Months with no clear scene are missing, so trend::sens.slope(post$nbr) (slope per STEP) would
  # mis-state the rate. The slope is the median of all pairwise slopes against real time (months
  # since the fire); the p-value is the Mann-Kendall test from trend::sens.slope (rank based).
  fire_date <- as.Date("2026-03-01")   # <-- REPLACE with the FIRMS date from Part A
  post <- dplyr::filter(recovery, date > fire_date) %>%
    mutate(months_since = as.numeric(date - fire_date) / 30.4375)
  ts_line <- NULL
  if (nrow(post) >= 4) {
    pw       <- utils::combn(nrow(post), 2)
    ts_slope <- median((post$nbr[pw[2, ]] - post$nbr[pw[1, ]]) /
                       (post$months_since[pw[2, ]] - post$months_since[pw[1, ]]))   # NBR per month
    ts_int   <- median(post$nbr - ts_slope * post$months_since)
    ts_p     <- trend::sens.slope(post$nbr)$p.value
    message("Post-fire NBR recovery (Theil-Sen): ", signif(ts_slope, 3), " NBR/month (",
            signif(ts_slope * 12, 3), " per year), p = ", signif(ts_p, 3), ", n = ", nrow(post))
    ts_line <- tibble(months_since = range(post$months_since)) %>%
      mutate(date = fire_date + months_since * 30.4375, nbr = ts_int + ts_slope * months_since)
    readr::write_csv(tibble(fire_date = fire_date, n_months = nrow(post), slope_per_month = ts_slope,
                            slope_per_year = ts_slope * 12, p_value = ts_p),
                     file.path(dir_tab, "nbr_recovery_trend.csv"))
  }

  rec_plot <- ggplot(recovery, aes(date, nbr)) +
    geom_vline(xintercept = fire_date, linetype = "dashed", colour = "red") +
    geom_line(colour = "grey60") + geom_point(colour = scico::scico(1, palette = "bamako")) +
    geom_line(data = ts_line, colour = "forestgreen", linewidth = 1) +   # no-op if ts_line is NULL
    labs(title = "Kitobo Forest — post-fire NBR recovery",
         subtitle = "Dashed line = fire date; green = Theil-Sen post-fire recovery trend",
         x = NULL, y = "Mean NBR (forest)") +
    theme_minimal(base_size = 12)
  print(rec_plot)
  ggsave(file.path(dir_fig, "03_regeneration_nbr.png"), rec_plot,
         width = 8, height = 4.5, dpi = 300)
  readr::write_csv(recovery, file.path(dir_tab, "nbr_recovery_series.csv"))
} else {
  message("Part C skipped: ", nbr_stack_tif, " not found. Run scripts/gee/02 and download to data/raw/gee/.")
}
