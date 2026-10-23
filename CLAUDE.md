# CLAUDE.md — Kitobo Forest Conservation Plan

> Project instructions for Claude Code. Read this fully before writing or running code.
> This file captures decisions already made with the project lead (Isaac Mureithi) so you
> can continue the work rather than re-deriving it. Last updated 2026-10-23, after the first
> hands-on pass: §2 records two corrections to the first handoff and §6 records what is verified.

---

## 1. What this project is

A spatial-analysis pipeline in **R (+ Google Earth Engine)** supporting the restoration and
post-fire monitoring of **Kitobo Forest** — a small groundwater/riverine forest in the Tsavo
landscape near Taveta, Taita-Taveta County, Kenya (AWF / OPALS work).

**Two analytical objectives:**

1. **Surface hydrology** — map "river starting points" (channel initiation / heads) and flow
   paths across the forest and its contributing catchment, plus a topographic wetness surface.
2. **Fire & regeneration** — map the extent and severity of the 2026 forest fire, and set up a
   vegetation-recovery time series to monitor regeneration and long-term restoration impact.

The repo is `Ismurn/Mureithi_AWF-Kitobo-Forest-Conservation-Plan`. It was created from the TESS Lab
`research-repository-template`, so the README follows that template and the layout is the lab standard
(`scripts/`, `data/raw`, `data/processed`, `outputs/`, `docs/`; GPL-3 licence).

---

## 2. Locked decisions (do NOT silently change these)

These were reasoned through already. If you believe one is wrong, flag it in your response
and explain — don't just override it.

| # | Decision | Why |
|---|----------|-----|
| D1 | **Survey points are in Arc 1960 / UTM Zone 37S (`EPSG:21037`); the working CRS for all metric work (area, buffers, hydrology, GEE) is WGS84 / UTM Zone 37S (`EPSG:32737`).** Set by `crs_pts` and `crs_utm` in `scripts/00_setup.R`. *(Revised 2026-10-23, see below.)* | The coordinates are projected UTM (eastings ~345 km, northings ~9,619 km = Taveta), NOT geographic. Other project layers are `EPSG:4326`; transform to match as needed. |
| D2 | **Build the forest polygon with a CONCAVE HULL (`sf::st_concave_hull`) of the survey points, NOT by joining the CSV points in file order.** Default `params$boundary = "hull"`, `params$hull_ratio = 0.10`. *(Reconfirmed by the lead on 2026-10-23, after seeing the comparison with the traced outline; see below.)* | The 74 points (70 unique) are not an ordered perimeter trace; joining in file order self-intersects (invalid geometry) and undercounts area at ~96 ha. Convex hull overshoots at ~200 ha. Concave hull `ratio≈0.1` ≈ 155 ha, `ratio≈0.3` ≈ 164 ha. Keep `ratio` a top-of-script parameter. Kitobo is usually cited ~150 ha. |
| D3 | **Hydrology = DEM + WhiteboxTools** (breach → D8 pointer → flow accumulation → extract streams → Strahler order → vectorise; TWI from SCA + slope). Slope/aspect/TPI/TRI via `terra::terrain`. | Surface channels are not visible from satellite here. This reuses the pipeline the lead already built for the Tsavo West rhino work. |
| D4 | **Run hydrology on a BUFFERED catchment (default 5 km), not the forest polygon alone.** | The forest sits at the discharge end; flow paths originate upslope. 155 ha of flat ground gives the DEM almost nothing to route. |
| D5 | **Fire date → NASA FIRMS; fire extent & severity → Sentinel-2 dNBR (GEE).** MODIS/VIIRS burned-area products (500 m) are too coarse — do not use them for extent. | FIRMS VIIRS is 375 m (1–3 pixels over this forest): good for *when*, useless for *where*. Sentinel-2 (10–20 m) dNBR gives the burn polygon + severity classes (Key & Benson). |
| D6 | **Regeneration = post-fire NBR/NDVI recovery time series**, trend via `trend::sens.slope()` (Theil-Sen), the lead's established method. | Defensible "regeneration rate" number for reporting. |

### Revisions after the first handoff (2026-10-23)

**D1 revised.** The first handoff assumed the CSV was WGS84 / UTM 37S. It is Arc 1960 / UTM 37S. Evidence:
against the operator's Google Earth outline (`Final Kitobo Forest.kml`, WGS84 by definition) the 70 unique points
have a median distance of **2.7 m** (90th percentile 3.5 m) when read as `EPSG:21037`, but **249 m** (317 m) when read
as `EPSG:32737`. The gap equals the Arc 1960 → WGS84 datum shift here (+97 m E, −299 m N): read as WGS84 the forest is
~315 m too far NNW. `01_boundary.R` repeats this check on every run and warns if another CRS fits better. To revert:
set `crs_pts <- "EPSG:32737"`. Figure: `docs/figures/kitobo_crs_check.png`.

**D2 reconfirmed by the lead (2026-10-23).** The same KML is an *ordered* trace of the perimeter. Closed, noded and
polygonized it gives one valid polygon of **176.5 ha** (`params$boundary = "outline"`). Against it, the hull at the
starting ratio 0.10 is 154.6 ha with IoU 0.79 and misses 30 ha of the outline; the best hull is ratio 0.50
(187.7 ha, IoU 0.92, but 12.9 ha beyond the outline). The ~150 ha figure matches a hull at ratio 0.05–0.10, not the
traced outline. The lead saw this comparison and chose the **concave hull**; the ratio stays 0.10 (it reproduces the ~150 ha usually quoted) unless they say otherwise, and the outline stays available for comparison. Say in reports that the hull covers about 83% of the traced outline. Table: `outputs/tables/01_hull_ratio_sensitivity.csv`;
figure: `docs/figures/kitobo_hull_vs_outline.png`. The KML is now only a QC reference (the boundary can be rebuilt from the CSV alone), so it stays gitignored; if the outline were ever chosen, `data/raw/Final_Kitobo_Forest.kml` would need to
be committed (it is gitignored for now) or the polygon saved, otherwise the repo cannot rebuild the boundary.

**Honest caveats to preserve in any report/README (these are scientific integrity points, not disclaimers to drop):**
- Modelled streams are **surface runoff/erosion drainage**, NOT the springs that actually feed
  this groundwater forest. Frame them as such; foreground **TWI** for the water-availability story,
  and pair with any field-known spring locations.
- A 30 m-class DEM over flat terrain makes channel-head location sensitive to the stream-extraction
  threshold and to canopy bias. Tune the threshold against imagery. For publication, prefer **FABDEM**
  (Copernicus GLO-30 with canopy/building bias removed) over `elevatr`'s blended tiles.
- dNBR quality depends on a reasonably cloud-free Sentinel-2 pre/post pair bracketing the FIRMS date.
  If the post window is cloudy, widen it or step to the next clear scene.

---

## 3. Tech stack & conventions (match the lead's existing style)

- **Package loading:** `pacman::p_load(...)`. Core spatial stack: `sf`, `terra`, `tidyverse`, `tmap` (v4 syntax — `tm_scale_*`, `col.scale`, `tm_title`; basemaps via `tm_basemap(server = ...)`, **not** `tm_view(basemap.server = ...)`, which v4 does not have). Zonal stats: `exactextractr`. Viz: `scico`, `viridis`, `patchwork`. Trend: `trend::sens.slope`. Hydrology: `whitebox`. DEM: `elevatr` (and FABDEM for publication).
- **Object naming:** dot notation — `Kitobo.Forest`, `Kitobo.DEM`, `Kitobo.Streams`, `Kitobo.TWI`.
- **Code style:** heavily commented, organised in clearly numbered **phases** (`# ---- 3b. ... ----`).
- **Paths:** relative to the repo root (open `research-repository-template.Rproj`, or run `Rscript` from the root). Machine-specific roots live only in the config block of `scripts/00_setup.R`.
- **GEE pattern:** author in the Code Editor → `Export.image.toDrive` → download the GeoTIFF into `data/raw/gee/` → import to R with `terra::rast()`. This mirrors the existing PML_V2 workflow. GEE account is **`ee-imnduta1`** (legacy `users/imnduta1`); scripts are saved under `users/imnduta1/Ismurn`.
- **WhiteboxTools is file-based:** write inputs to disk (`data/processed/wbt/`), pass file paths to `wbt_*`, read results back. The binary is installed (v2.4.0); `wbt_init()` runs in `00_setup.R`.
- **Machine notes (lead's PC):** `Rscript` is not on PATH — use `C:\Program Files\R\R-4.4.1\bin\Rscript.exe`. Pandoc exists only inside RStudio, so `save_map()` writes non-self-contained HTML (an `_files` folder beside it) from plain `Rscript`. To knit the results note from plain `Rscript` (as tested), first `Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")` and pass `output_dir` to `rmarkdown::render()`, which otherwise writes next to the Rmd. The session scratch path is longer than Windows' 260-character limit, so put temporary scripts in `%TEMP%`, not in the session folder.

---

## 4. Data inventory

| Layer | Location | Notes |
|-------|----------|-------|
| Kitobo survey points | `data/raw/Kitobo_Forest_Coordinates.csv` | 74 rows (70 unique), columns `X,Y` in **`EPSG:21037`** (Arc 1960 / UTM 37S). Same as the lead's `Kitobo Forest  Coordinates.csv` apart from line endings. |
| Google Earth outline | `data/raw/Final_Kitobo_Forest.kml` (gitignored; original in the lead's OneDrive, `R WORK/CommunityForestAssociations-Tsavo/Kitobo Forest/`) | One LineString, 102 vertices, WGS84. CRS reference and a comparison for the hull (D2); not used for the boundary. |
| Kitobo Forest polygon | `outputs/shapefiles/Kitobo_Forest.shp` (+ `data/processed/kitobo_forest.gpkg`) | **Produced by `scripts/01_boundary.R`**, EPSG:32737. Upload to GEE as `projects/ee-imnduta1/assets/Kitobo_Forest` for the fire scripts. |
| Tsavo ranches / settled areas | `C:/Workspace/Mureithi_RPI_Conservation_Governance/data/Tsavo_Conservation_Areas` (`dir_tsv`) | `RanchesTTWCA.shp`, `SettledAreas.shp`. Copies also sit in two other repos under `C:/Workspace`. Kitobo lies **outside** both layers (nearest feature ~55 km east), so there is no host ranch to report. |
| Kenya admin / protected areas | `C:/Workspace/Mureithi_RPI_Conservation_Governance/data/Country_Borders` (`dir_kenya`) | `*_wgs84.shp` (EPSG:4326). The admin layer is the national outline only. Nearest protected area: Tsavo West NP, 19.2 km (`outputs/tables/01_landscape_context.csv`). |
| GEE dNBR / severity / NBR rasters | `data/raw/gee/` | `Kitobo_dNBR.tif` and `Kitobo_Severity.tif` (10 m, 624 × 511), `Kitobo_NBR_monthly.tif` (20 m, 13 monthly bands Oct 2025–Oct 2026), all float32 in EPSG:32737 on the hull's 2 km buffer, from the lead's run of 2026-10-23 (the Drive zip is kept beside them). Imported in `scripts/03_fire_regeneration.R`. |

External data lives on the lead's machine (per existing convention); do not assume absolute paths —
read roots from the config block at the top of `scripts/00_setup.R`.

---

## 5. Repo layout

```
.
├── CLAUDE.md                     # this file
├── README.md                     # project overview (TESS Lab README template, filled in)
├── LICENSE                       # GPL-3 (from the template)
├── .gitignore                    # template rules + secrets, data/processed, outputs, rasters
├── research-repository-template.Rproj
├── scripts/
│   ├── 00_setup.R                # packages, CRS, paths/config, params, helpers
│   ├── 01_boundary.R             # CRS check → boundary polygon → landscape context → shapefile
│   ├── 02_dem_hydrology.R        # DEM over buffered catchment; terrain; flow/streams/heads/TWI
│   ├── 03_fire_regeneration.R    # FIRMS query; import dNBR/severity; burn stats; recovery + sens.slope
│   ├── gee/
│   │   ├── 01_kitobo_dNBR.js         # Sentinel-2 burn severity (dNBR) → Drive
│   │   ├── 02_kitobo_regeneration.js # NBR/NDVI recovery chart + monthly NBR stack → Drive
│   │   └── kitobo_gee_combined_as_run_2026-10-23.js  # the lead's combined script as run; produced data/raw/gee/
│   └── Mureithi_AWF-Kitobo-Forest-Conservation-Plan.Rmd  # results note (knit to HTML or Word); reads the pipeline outputs
├── data/
│   ├── raw/                      # Kitobo_Forest_Coordinates.csv, Final_Kitobo_Forest.kml, gee/ (GEE downloads)
│   └── processed/                # intermediates written by the pipeline, wbt/ scratch
├── outputs/
│   ├── shapefiles/  figures/  rasters/  tables/
└── docs/
    ├── methodology.md            # reasoned write-up + caveats (keep in sync with decisions)
    └── figures/                  # figures referenced by methodology.md
```

---

## 6. Build roadmap — status and known issues

Status as of 2026-10-23. `[x]` done and verified, `[~]` runs but needs a decision or tuning, `[ ]` not started.

- [x] **T1 — Environment.** R 4.4.1; `sf` 1.0-19 (GEOS 3.12.2, so `st_concave_hull` works), `terra` 1.8-54, `tmap` 4.1, `whitebox` 2.4.3 with WhiteboxTools v2.4.0 already installed, `elevatr` 0.99.1, `trend` 1.1.6. `exactextractr` 0.10.1 was missing and was installed from CRAN (`pacman::p_load` does this on first run). **Not done:** `renv` — no lockfile yet; the lead's call (the TESS template recommends it, and initialising it changes how R starts in this project).
- [~] **T2 — Boundary (`01_boundary.R`).** Runs end to end in ~25 s. Writes `Kitobo_Forest.shp` (hull, 154.6 ha by default), the CRS check, the hull sensitivity table, two comparison figures and the landscape context. **Decided:** D2 stays the hull at ratio 0.10 (lead, 2026-10-23). **Open:** an imagery check by the lead (`outputs/figures/01_boundary_qc.html`). **Not possible as written:** "confirm the host ranch / settlement" — Kitobo is outside the Tsavo layers (§4).
- [~] **T3 — Hydrology (`02_dem_hydrology.R`).** Runs end to end in ~3 min (needs internet for the DEM tiles). The first run showed that the scaffold's 1,000-cell stream threshold was only 0.09 km² on the 9.5 m grid `elevatr` returns (interpolated from ~30 m data): 613 segments, 321 km, 2.9 km/km², 327 channel heads — mostly DEM noise on a flat plain. The threshold is now `params$stream_area_km2`, converted to cells with the real resolution, provisionally **1 km²**: 58 segments, 98 km, 0.9 km/km², 35 channel heads, none inside the forest (`outputs/figures/02_hydrology.png` and `.html`; TWI in `outputs/rasters/kitobo_twi.tif`). Checks passed: streams, DEM and TWI share CRS (32737) and extent. **Open:** tune `stream_area_km2` against imagery and any field-known channels, and decide whether to move to FABDEM before publication.
- [ ] **T4 — Fire date (`03` Part A).** *Parked by the lead on 2026-10-23: the FIRMS check is not needed for now, so the fire date stays "mid-February 2026" as given.* If it is ever wanted: wire up `FIRMS_MAP_KEY` (§7), run, and check the fire date against the lead's "mid-February 2026" (`params$fire_date` = 2026-02-15; it only needs the right month because the composites are monthly, and the February composite already shows the burn). **Fixed but untested with a valid key:** the first scaffold used DAY_RANGE 10, but the API allows **1–5**, and its `try()` would have turned every failure into "no fires". Part A now sweeps in 5-day steps, queries `_SP` and `_NRT` VIIRS sources, de-duplicates, and reports error replies (a bad key gives `HTTP 400: Invalid MAP_KEY.`, checked). Once it runs, consider counting detections within ~1 km of the forest (the current AOI is the 5 km buffer) to separate the forest fire from regional burning, and update `params$fire_date` if FIRMS disagrees.
- [x] **T5 — dNBR (Earth Engine → `03` Part B).** Done 2026-10-23 on the lead's exports (`data/raw/gee/GEE_Exports-20261023T122407Z-1-001.zip`, unpacked beside it). The script they ran is `scripts/gee/kitobo_gee_combined_as_run_2026-10-23.js` (a combined 01 + 02; windows 10 Jan–10 Feb vs 1 Mar–15 Apr 2026; asset `users/imnduta1/Kitobo_Forest`, which matches the hull's 2 km buffer to one pixel). The exports are float32 GeoTIFFs in EPSG:32737 (dNBR and severity 10 m, 624 × 511; monthly NBR 20 m, 312 × 256, 13 bands named `YYYY_MM_dd_NBR`, so the band-name fix worked). **Result:** 23.5% of the hull (36.3 ha) is burned (dNBR > 100): Low 12.9%, Moderate-low 9.0%, Moderate-high 1.3%, High 0.3% of the forest, so 93% of the burned area is low or moderate-low. Under the traced outline it is 35.4 ha (20.1%): the hectares barely depend on the boundary, only the share does (`outputs/tables/burn_by_boundary.csv`; the hull's 30 ha "missed" wedge holds 0.3 ha of burn). Inside the forest the unburned ground is steady (median dNBR −2), so the Key & Benson classes are usable there; outside it the savanna greens up between the windows (median −109), so never read severity outside the forest. Part B now rounds the float severity to whole classes.
- [x] **T6 — Regeneration (Earth Engine → `03` Part C).** Done 2026-10-23 on the real monthly NBR stack (all 13 months have data over the forest). The scaffold's whole-forest series was the wrong measure: only 23.5% of the forest burned and the series carries the seasonal cycle of the rains, so its Theil-Sen slope after the fire is flat (−0.0006 NBR/month, p = 0.90). Part C now measures the **gap between burned and unburned ground inside the forest** (same vegetation and season; burned = at least half of the four 10 m cells in a 20 m cell burned, unburned = none, mixed cells left out). The gap was +0.002 before the fire (so the unburned ground is a fair control), fell to −0.36 in March 2026 and was within 10% of its pre-fire value by June 2026, 3.5 months after the fire. The Theil-Sen slope over that recovery phase (March–June, n = 4) is +0.125 NBR/month (p = 0.09: four points give a test with little power). Slopes use real time (months since `params$fire_date`), not per-step. Outputs: `nbr_recovery_series.csv`, `nbr_recovery_trend.csv`, `nbr_recovery_summary.csv`, `03_regeneration_nbr.png`. **Caveats to keep:** this is spectral recovery, in step with the long rains, not woody structure, species or springs; the series has only four pre-fire months, so seasonality-adjusted trends need a multi-year baseline. **To do:** extend the export as new months arrive; add field indicators of structural recovery.
- [~] **T7 — Reporting.** `scripts/Mureithi_AWF-Kitobo-Forest-Conservation-Plan.Rmd` (the lead's empty notebook) is now the results note: key-results table, boundary and CRS check, hull comparison, hydrology and TWI, fire date, burn severity and burned area by boundary, regeneration (burned vs unburned), caveats, package versions. It reads `outputs/` (tables, figures, rasters) instead of recomputing the pipeline, so run the numbered scripts first; sections without data print "Pending". Knit to HTML and Word with the real data on 2026-10-23 (`outputs/Mureithi_AWF-Kitobo-Forest-Conservation-Plan.html` and `.docx`, gitignored). Chunks can also be run one by one in RStudio, which starts in `scripts/`: the first chunk sets the working directory to the repo root when not knitting and `setup` sources `here::here("scripts", "00_setup.R")` (the first version failed with "cannot open file 'scripts/00_setup.R'"; tested by running all chunks from `scripts/` and from the root). **To do:** re-knit after any rerun. **Parked by the lead (2026-10-23):** matching the AWF/donor table layout (they do not have it yet; a `reference_docx` template would carry the Word styling) and the FIRMS check of the fire date (T4).

---

## 7. Secrets & what not to commit

- **FIRMS `MAP_KEY`** (free, from https://firms.modaps.eosdis.nasa.gov/api/) goes in a **`.Renviron`**
  file at repo root as `FIRMS_MAP_KEY=...`. `.Renviron` is gitignored (verified). Never hardcode the key in a
  script or commit it, and never ask for it in chat — the lead creates the file. A plain `Rscript` started
  outside the repo root will not read it; call `readRenviron(".Renviron")` first.
- `data/processed/`, `outputs/`, `data/raw/gee/` and all `*.tif` are gitignored (large/derived; `dir/*` rules so
  the `.gitkeep` files still work). `data/raw/Kitobo_Forest_Coordinates.csv` is tracked (the lead committed
  and pushed it on 2026-10-23, commit `7b0b94d`, which did not include `.gitignore` or `CLAUDE.md`). `data/raw/*.kml` is
  ignored on purpose.
- Do not commit anything from the lead's OneDrive data roots.

---

## 8. Reproduce / run order

```bash
# from the repo root, in R (or Rscript), in this order — each script sources 00_setup.R itself:
Rscript scripts/01_boundary.R          # CRS check, boundary polygon, landscape context, shapefile
Rscript scripts/02_dem_hydrology.R     # needs internet for DEM tiles
# --- upload outputs/shapefiles/Kitobo_Forest.shp to GEE asset projects/ee-imnduta1/assets/Kitobo_Forest ---
Rscript scripts/03_fire_regeneration.R # first pass: Part A only (FIRMS date; needs FIRMS_MAP_KEY)
# --- run the Earth Engine script (scripts/gee/; the lead's as-run combined version) with windows bracketing the fire; download exports to data/raw/gee/ ---
Rscript scripts/03_fire_regeneration.R # second pass: Parts B and C
```

Scripts write intermediates to `data/processed/` and `outputs/`, so later scripts reload what they need;
but running in one RStudio session also works and is what the lead usually does.
