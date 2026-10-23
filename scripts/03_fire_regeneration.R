# 03_fire_regeneration.R  —  Fire extent/severity + regeneration monitoring
# Part A: FIRMS sweep to find the 2026 fire date(s).
# Part B: import Sentinel-2 dNBR/severity (from scripts/gee/01_kitobo_dNBR.js) and map.
# Part C: regeneration — NBR/NDVI recovery series + Theil-Sen trend.
# Prereq: 01_boundary.R (polygon). Parts B/C need the GEE exports in data/raw/gee/.
source("scripts/00_setup.R")

Kitobo.Forest <- st_read(file.path(dir_proc, "kitobo_forest.gpkg"), quiet = TRUE) %>%
  st_transform(crs_utm)
Kitobo.Buff <- st_buffer(Kitobo.Forest, dist = params$catchment_buffer)

# PART A — FIRMS: find WHEN it burned (375 m VIIRS: good for date, not extent)
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

# PART B — Sentinel-2 dNBR: WHERE & HOW SEVERE (import GEE export, then map)
# Run scripts/gee/01_kitobo_dNBR.js first; download the GeoTIFFs into data/raw/gee/.
dnbr_tif <- file.path(dir_gee, "Kitobo_dNBR.tif")
sev_tif  <- file.path(dir_gee, "Kitobo_Severity.tif")

if (file.exists(sev_tif)) {
  Kitobo.dNBR <- rast(dnbr_tif)
  Kitobo.Sev  <- round(rast(sev_tif))   # the GEE export stores the classes as float32; round to whole classes

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

  # Burned area (dNBR > 100) under each boundary on offer: the forest polygon in use and, if 01_boundary.R
  # saved it, the operator's traced outline. The hectares barely move with the boundary; the share does.
  zones <- list("Forest boundary used" = Kitobo.Forest)
  if (file.exists(file.path(dir_proc, "kitobo_outline.gpkg")))
    zones[["Traced Google Earth outline"]] <- st_read(file.path(dir_proc, "kitobo_outline.gpkg"), quiet = TRUE)
  burned_ras <- (Kitobo.Sev >= 1) * 1
  burn_by_boundary <- purrr::imap_dfr(zones, function(z, nm) {
    z     <- st_transform(z, crs(Kitobo.Sev))
    share <- as.numeric(exact_extract(burned_ras, z, fun = "mean"))
    tibble(boundary = nm, area_ha = area_ha(z), burned_share_pct = 100 * share, burned_ha = area_ha(z) * share)
  })
  readr::write_csv(burn_by_boundary, file.path(dir_tab, "burn_by_boundary.csv"))
  print(burn_by_boundary)
} else {
  message("Part B skipped: ", sev_tif, " not found. Run the Earth Engine script (scripts/gee/) and download the exports to data/raw/gee/.")
}

# PART C — Regeneration: monthly NBR of burned vs unburned forest + Theil-Sen recovery rate
# Needs the monthly NBR export (scripts/gee/) in data/raw/gee/.
# Whole-forest NBR mixes burned and unburned ground with the seasonal swing of the rains, so the
# recovery is measured as the GAP between burned and unburned pixels INSIDE the forest (same
# vegetation, same season): near 0 before the fire, a drop at the fire, back to ~0 as the burned
# ground recovers. Burned = dNBR > 100 (Part B severity classes 1 to 4).
nbr_stack_tif <- file.path(dir_gee, "Kitobo_NBR_monthly.tif")   # bands = months, named YYYY_MM_dd_NBR

if (file.exists(nbr_stack_tif)) {
  NBR.Stack <- rast(nbr_stack_tif)
  NBR.Stack[NBR.Stack == 0] <- NA   # GEE can write masked (cloudy) pixels as 0; real NBR is never exactly 0

  # Mean NBR over the whole forest per month (exact_extract weights each cell by the area inside the polygon)
  forest_nbr <- st_transform(Kitobo.Forest, crs(NBR.Stack))
  mean_nbr   <- function(stack) as.numeric(exact_extract(stack, forest_nbr, fun = "mean"))
  dates <- as.Date(substr(names(NBR.Stack), 1, 10), format = "%Y_%m_%d")   # band names start with the month
  if (all(is.na(dates)))
    stop("Band names are not dates: re-export with the monthly NBR script (names must start YYYY_MM_dd).")
  recovery <- tibble(date = dates, nbr_forest = mean_nbr(NBR.Stack))

  # Burned vs unburned ground INSIDE the forest, from the 10 m severity raster resampled to the 20 m NBR
  # grid: burned = at least half of the 4 sub-cells burned; unburned = none of them; mixed cells are left
  # out of both. Before the fire the two should sit together (the check that unburned is a fair control).
  if (file.exists(sev_tif)) {
    burn_share <- resample((round(rast(sev_tif)) >= 1) * 1, NBR.Stack[[1]], method = "average")
    recovery <- recovery %>% mutate(
      nbr_burned   = mean_nbr(mask(NBR.Stack, ifel(burn_share >= 0.5, 1, NA))),
      nbr_unburned = mean_nbr(mask(NBR.Stack, ifel(burn_share == 0, 1, NA))),
      gap          = nbr_burned - nbr_unburned)
  }
  recovery <- arrange(recovery, date)

  # Theil-Sen rate of change against REAL time (months since the fire). Months without a clear scene are
  # missing, so trend::sens.slope() on the bare series (slope per step) would mis-state the rate. The slope
  # is the median of all pairwise slopes; the p-value is the Mann-Kendall test from trend::sens.slope().
  fire_date    <- params$fire_date                                  # set in 00_setup.R; refine from Part A
  months_since <- function(d) as.numeric(d - fire_date) / 30.4375
  theil_sen <- function(measure, d) {                               # d = tibble(date, y)
    t <- months_since(d$date)
    s <- median(utils::combn(nrow(d), 2, function(i) (d$y[i[2]] - d$y[i[1]]) / (t[i[2]] - t[i[1]])))
    tibble(measure = measure, from = min(d$date), to = max(d$date), n_months = nrow(d),
           slope_per_month = s, slope_per_year = 12 * s,
           p_value = if (nrow(d) >= 4) trend::sens.slope(d$y)$p.value else NA_real_,
           intercept = median(d$y - s * t))
  }

  trends   <- tibble()
  rec_sum  <- NULL
  post_all <- recovery %>% filter(date > fire_date, !is.na(nbr_forest)) %>% transmute(date, y = nbr_forest)
  if (nrow(post_all) >= 4)
    trends <- bind_rows(trends, theil_sen("Whole-forest NBR, all post-fire months", post_all))

  if ("gap" %in% names(recovery) && sum(!is.na(recovery$gap[recovery$date > fire_date])) >= 2) {
    pre     <- filter(recovery, date < floor_date(fire_date, "month"), !is.na(gap))   # months wholly before the fire month
    post    <- filter(recovery, date > fire_date, !is.na(gap))
    gap_pre <- mean(pre$gap)
    trough  <- post[which.min(post$gap), ]
    thr     <- gap_pre + 0.10 * (trough$gap - gap_pre)                                # within 10% of the pre-fire gap = 90% recovered
    recov   <- post %>% filter(date >= trough$date, gap >= thr) %>% slice_head(n = 1)
    phase   <- post %>% filter(date >= trough$date, date <= if (nrow(recov)) recov$date else max(date))
    if (nrow(phase) >= 2)
      trends <- bind_rows(trends, theil_sen("Gap burned minus unburned, trough to recovery", transmute(phase, date, y = gap)))
    rec_sum <- tibble(fire_date = fire_date, pre_fire_gap = gap_pre, trough_month = trough$date, trough_gap = trough$gap,
                      recovered_month = if (nrow(recov)) recov$date else as.Date(NA),
                      months_fire_to_recovered = if (nrow(recov)) months_since(recov$date) else NA_real_)
    message("Burned-minus-unburned NBR gap: ", round(gap_pre, 3), " before the fire, trough ", round(trough$gap, 3),
            " in ", format(trough$date, "%b %Y"), "; ",
            if (nrow(recov)) paste0("within 10% of the pre-fire gap by ", format(recov$date, "%b %Y"), " (",
                                    round(months_since(recov$date), 1), " months after the fire)") else "not yet recovered")
  }
  for (i in seq_len(nrow(trends)))
    message("Theil-Sen, ", trends$measure[i], ": ", signif(trends$slope_per_month[i], 3), " NBR/month (",
            signif(trends$slope_per_year[i], 3), " per year), n = ", trends$n_months[i],
            ", p = ", signif(trends$p_value[i], 3))

  # Figure: monthly NBR of the three series, then (if severity is available) the gap with its recovery slope
  long <- recovery %>%
    select(date, any_of(c("nbr_forest", "nbr_burned", "nbr_unburned"))) %>%
    pivot_longer(-date, names_to = "series", values_to = "nbr") %>%
    mutate(series = recode(series, nbr_forest = "Whole forest", nbr_burned = "Burned ground",
                           nbr_unburned = "Unburned ground"))
  fire_line <- geom_vline(xintercept = fire_date, linetype = "dashed", colour = "red")
  rec_plot  <- ggplot(long, aes(date, nbr, colour = series)) + fire_line +
    geom_line() + geom_point(size = 1.6) +
    scale_colour_manual(NULL, values = c("Whole forest" = "grey50", "Burned ground" = "firebrick",
                                         "Unburned ground" = "forestgreen")) +
    labs(title = "Kitobo Forest: monthly NBR after the 2026 fire",
         subtitle = "Dashed line = fire date (mid-month)", x = NULL, y = "Mean NBR") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  if ("gap" %in% names(recovery)) {
    gap_plot <- ggplot(recovery, aes(date, gap)) + fire_line + geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(colour = "grey40") + geom_point(colour = "firebrick", size = 1.6) +
      labs(subtitle = "Gap = burned minus unburned NBR; green = Theil-Sen slope from the trough to recovery",
           x = NULL, y = "NBR gap") +
      theme_minimal(base_size = 11)
    g <- trends[grepl("^Gap", trends$measure), ]
    if (nrow(g))
      gap_plot <- gap_plot + geom_line(
        data = tibble(date = c(g$from, g$to)) %>% mutate(gap = g$intercept + g$slope_per_month * months_since(date)),
        colour = "forestgreen", linewidth = 1.1)
    rec_plot <- rec_plot / gap_plot
  }
  if (interactive()) print(rec_plot)   # (printing under Rscript would leave an Rplots.pdf behind)
  ggsave(file.path(dir_fig, "03_regeneration_nbr.png"), rec_plot, width = 8,
         height = if ("gap" %in% names(recovery)) 7 else 4.5, dpi = 300)
  readr::write_csv(recovery, file.path(dir_tab, "nbr_recovery_series.csv"))
  if (nrow(trends)) readr::write_csv(select(trends, -intercept), file.path(dir_tab, "nbr_recovery_trend.csv"))
  if (!is.null(rec_sum)) readr::write_csv(rec_sum, file.path(dir_tab, "nbr_recovery_summary.csv"))
} else {
  message("Part C skipped: ", nbr_stack_tif, " not found. Run the Earth Engine script (scripts/gee/) and download the exports to data/raw/gee/.")
}
