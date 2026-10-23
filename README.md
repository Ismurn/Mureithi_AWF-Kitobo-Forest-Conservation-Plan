# Kitobo Forest Conservation Plan

This repository contains the data and code for the spatial analysis supporting the restoration and post-fire monitoring of **Kitobo Forest**, a small groundwater/riverine forest in the Tsavo landscape near Taveta, Taita-Taveta County, Kenya (African Wildlife Foundation / OPALS).

Use of this code is licensed under [![Licence: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Contact Isaac Mureithi at *IMureithi@awf.org / I.N.Mureithi@exeter.ac.uk*. Citation and Zenodo DOI: to be added on publication.

---

## Project overview

The analysis has two objectives:

1. **Surface hydrology**: map river/channel starting points and flow paths across the forest and its contributing catchment, with a topographic wetness surface (TWI).
2. **Fire and regeneration**: map the extent and severity of the 2026 forest fire and track vegetation recovery to assess long-term restoration impact.

It is an R pipeline (`sf`, `terra`, WhiteboxTools) with Google Earth Engine scripts for the Sentinel-2 steps. The reasoning behind each method, and the caveats that must travel with the results, are in [*docs/methodology.md*](docs/methodology.md). Two points matter most when reading the outputs:

- Modelled streams are **surface runoff/erosion drainage**, not the springs that feed this groundwater forest. The TWI surface is the better guide to water availability.
- The surveyed boundary coordinates are in **Arc 1960 / UTM 37S** (EPSG:21037), not WGS84. Read as WGS84 the forest is displaced by about 315 m.

---

## Repository contents

### Scripts

| Script | Description |
|---------|-------------|
| *scripts/00_setup.R* | Packages, CRS, paths, tunable parameters (`params`) and helpers. Sourced by every phase script. |
| *scripts/01_boundary.R* | Reads the survey points, checks their CRS against the Google Earth outline, builds the forest polygon, reports landscape context and saves *Kitobo_Forest.shp*. |
| *scripts/02_dem_hydrology.R* | DEM over a buffered catchment, terrain derivatives, WhiteboxTools flow routing, stream network, channel heads and TWI. |
| *scripts/03_fire_regeneration.R* | Part A: NASA FIRMS sweep for the fire date. Part B: import Sentinel-2 dNBR severity and map it. Part C: NBR recovery series and Theil-Sen trend. |
| *scripts/gee/01_kitobo_dNBR.js* | Earth Engine: Sentinel-2 dNBR and Key & Benson severity classes, exported to Drive. |
| *scripts/gee/02_kitobo_regeneration.js* | Earth Engine: NBR/NDVI recovery chart and monthly NBR stack, exported to Drive. |
| *scripts/Mureithi_AWF-Kitobo-Forest-Conservation-Plan.Rmd* | Results note (knits to HTML or Word). It reads the pipeline outputs, so the boundary and hydrology sections fill in now and the fire and regeneration sections fill in once *scripts/03_fire_regeneration.R* has run on the Earth Engine exports. |

### Data

| File | Description |
|------|-------------|
| *data/raw/Kitobo_Forest_Coordinates.csv* | 74 surveyed boundary points (70 unique); `X`,`Y` in Arc 1960 / UTM 37S (EPSG:21037). |
| *data/raw/Final_Kitobo_Forest.kml* | The operator's Google Earth outline of the forest (WGS84). Used as a CRS reference and to compare the hull with the traced outline. Kept locally and not tracked (see *.gitignore*); the pipeline runs without it. |
| *data/raw/gee/* | GeoTIFFs downloaded from Earth Engine (Drive exports). Not tracked: regenerable. |
| *data/processed/* | Intermediate files written by the pipeline (GeoPackage copies, WhiteboxTools scratch). Not tracked. |
| *outputs/* | Shapefiles, rasters, figures and tables written by the scripts. Not tracked. |

The Tsavo ranch and settlement layers (*RanchesTTWCA.shp*, *SettledAreas.shp*) and the Kenya administrative and protected-area layers are read from the folders set in `dir_tsv` and `dir_kenya` in *scripts/00_setup.R*. They are held in other TESS Lab repositories and are not redistributed here. Kitobo lies west of the ranch and settlement layers, so they only supply distances (*outputs/tables/01_landscape_context.csv*).

---

## Getting started

Clone this repository and open *research-repository-template.Rproj* in RStudio. The scripts use paths relative to the repository root.

- **R packages** are installed on first run by `pacman::p_load()` in *scripts/00_setup.R*. The pipeline was developed on R 4.4.1 with `sf` 1.0-19 (GEOS 3.12; `st_concave_hull()` needs GEOS 3.11 or later), `terra` 1.8-54, `tmap` 4.1, `whitebox` 2.4.3, `elevatr` 0.99.1 and `exactextractr` 0.10.1. No *renv.lock* has been created yet.
- **WhiteboxTools**: run `whitebox::install_whitebox()` once.
- **External data roots**: edit `dir_tsv` and `dir_kenya` in *scripts/00_setup.R* to match your machine.
- **NASA FIRMS key** (free, <https://firms.modaps.eosdis.nasa.gov/api/>): add the line `FIRMS_MAP_KEY=...` to a `.Renviron` file in the repository root. It is gitignored; never commit it.
- **Earth Engine**: the scripts are written for account `ee-imnduta1`.

## Running the analysis

Each R script sources *scripts/00_setup.R* itself. Run them in this order:

1. *scripts/01_boundary.R*. Then upload *outputs/shapefiles/Kitobo_Forest.shp* to Earth Engine as `projects/ee-imnduta1/assets/Kitobo_Forest`.
2. *scripts/02_dem_hydrology.R*
3. *scripts/03_fire_regeneration.R*, once, for Part A only (it reports the fire date from FIRMS).
4. *scripts/gee/01_kitobo_dNBR.js* and *scripts/gee/02_kitobo_regeneration.js* in the Earth Engine Code Editor, with the pre/post windows set from the FIRMS date. Download the Drive exports into *data/raw/gee/*.
5. *scripts/03_fire_regeneration.R* again, for Parts B and C (burn severity, recovery trend).
6. Knit *scripts/Mureithi_AWF-Kitobo-Forest-Conservation-Plan.Rmd* (HTML or Word) to assemble the results note. It can be knitted at any stage; sections without data show as pending.

---

## Methodological notes and caveats

- **CRS.** Survey points are Arc 1960 / UTM 37S; all metric work uses WGS84 / UTM 37S (EPSG:32737). *01_boundary.R* re-checks the points against the Google Earth outline on every run.
- **Boundary.** The forest polygon is the concave hull (ratio 0.10, 154.6 ha) of the survey points. Because the points lie only on the perimeter, a low ratio cuts chords across the interior: the hull covers about 83% of the 176 ha enclosed by the operator's traced outline. See *docs/methodology.md* section 2.
- **Hydrology.** A 30 m-class DEM over flat terrain makes channel-head position sensitive to the stream threshold and to canopy bias. FABDEM is preferred for publication.
- **Fire.** FIRMS (375 m) gives the date only; extent and severity come from Sentinel-2 dNBR (10-20 m). dNBR needs a reasonably cloud-free pre/post pair bracketing the fire date.

Work plan and status: see [*CLAUDE.md*](CLAUDE.md) section 6.
