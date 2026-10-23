// =============================================================================
// AS RUN by the lead in the Earth Engine Code Editor on 2026-10-23.
// A combined version of 01_kitobo_dNBR.js and 02_kitobo_regeneration.js, with the windows set for a
// mid-February 2026 fire. It produced data/raw/gee/Kitobo_dNBR.tif, Kitobo_Severity.tif and
// Kitobo_NBR_monthly.tif (Drive folder GEE_Exports, downloaded as GEE_Exports-20261023T122407Z-1-001.zip).
// Everything below this header is the lead's script, saved verbatim as the record of the run.
// =============================================================================

// ============================================================================
// Kitobo Forest Conservation Plan — Comprehensive GEE Script
// Integrates:
//   1. Burn severity mapping & classification (dNBR & Key & Benson)
//   2. Monthly NBR/NDVI time-series stack for R-side Theil-Sen regeneration analysis
// ============================================================================

// ---- 0. Area of interest ---------------------------------------------------
var kitobo = ee.FeatureCollection('users/imnduta1/Kitobo_Forest');
var aoi    = kitobo.geometry().buffer(2000);   // 2 km buffer to capture burn edges
Map.centerObject(aoi, 14);

// ---- 1. Fire-bracketing windows for dNBR (Mid-Feb 2026 fire) ---------------
var preStart  = '2026-01-10', preEnd  = '2026-02-10';   // Clear window BEFORE fire
var postStart = '2026-03-01', postEnd = '2026-04-15';   // Clear window AFTER fire

// ---- 2. Sentinel-2 Preprocessing Functions (Masking & Indices) ------------
function maskS2(img) {
  var scl = img.select('SCL');
  // Drop: 3 (shadow), 8/9 (cloud med/high), 10 (cirrus), 11 (snow)
  var ok = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(ok).divide(10000)
            .copyProperties(img, ['system:time_start']);
}

function addIndices(img) {
  var nbr  = img.normalizedDifference(['B8', 'B12']).rename('NBR');
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return img.addBands([nbr, ndvi]);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(aoi)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
           .map(maskS2).map(addIndices);

// ---- 3. dNBR + Key & Benson Severity Classes ------------------------------
var preNBR  = s2.filterDate(preStart,  preEnd ).select('NBR').median();
var postNBR = s2.filterDate(postStart, postEnd).select('NBR').median();

var dNBR = preNBR.subtract(postNBR).multiply(1000).rename('dNBR').clip(aoi);

var severity = ee.Image(0)
  .where(dNBR.gt(100).and(dNBR.lte(270)), 1)    // Low
  .where(dNBR.gt(270).and(dNBR.lte(440)), 2)    // Moderate-low
  .where(dNBR.gt(440).and(dNBR.lte(660)), 3)    // Moderate-high
  .where(dNBR.gt(660),                        4)    // High
  .rename('severity').clip(aoi).toByte();

// ---- 4. Monthly Regeneration Stack (Baseline to Present) ------------------
var startDate = '2025-10-01';
var endDate   = '2026-10-23';

var s2Monthly = s2.filterDate(startDate, endDate);

var months = ee.List.sequence(0,
  ee.Date(endDate).difference(ee.Date(startDate), 'month').round().subtract(1));
var blank = ee.Image.constant([0, 0]).rename(['NBR', 'NDVI']).selfMask();

var monthly = ee.ImageCollection.fromImages(months.map(function(m) {
  var t0 = ee.Date(startDate).advance(m, 'month');
  var t1 = t0.advance(1, 'month');
  var col  = s2Monthly.filterDate(t0, t1).select(['NBR', 'NDVI']);
  var comp = ee.Image(ee.Algorithms.If(col.size().gt(0), col.median(), blank));
  return comp.set({'system:time_start': t0.millis(),
                   'system:index': t0.format('YYYY_MM_dd')});
}));

// ---- 5. On-Screen Visualization & Console Inspections ---------------------
var sevPalette = ['d9d9d9', 'ffffb2', 'fecc5c', 'fd8d3c', 'e31a1c'];
Map.addLayer(dNBR, {min: -200, max: 800, palette: ['006400','ffffff','ff0000']}, 'dNBR');
Map.addLayer(severity, {min: 0, max: 4, palette: sevPalette}, 'Severity Classes');
Map.addLayer(kitobo.style({color: '000000', fillColor: '00000000', width: 2}), {}, 'Kitobo Boundary');

// Recovery trajectory chart
var chart = ui.Chart.image.series({
  imageCollection: monthly.select(['NBR', 'NDVI']),
  region: kitobo.geometry(), reducer: ee.Reducer.mean(), scale: 20, xProperty: 'system:time_start'
}).setOptions({
  title: 'Kitobo Forest — Monthly Mean NBR & NDVI Recovery',
  hAxis: {title: 'Date'}, vAxis: {title: 'Index Value'}, lineWidth: 2, pointSize: 3
});
print(chart);

// Burned area breakdown printed to console
var areaImg = ee.Image.pixelArea().divide(1e4).addBands(severity);
var areas = areaImg.reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'severity'}),
  geometry: kitobo.geometry(), scale: 10, maxPixels: 1e9
});
print('Burned area (ha) by severity class within Kitobo:', areas);

// ---- 6. Google Drive Exports (Using robust export structure) --------------
// 1. dNBR Raster Export
Export.image.toDrive({
  image: dNBR.toFloat(),
  description: 'Kitobo_dNBR',
  folder: 'GEE_Exports',
  region: aoi.bounds(),
  scale: 10,
  crs: 'EPSG:32737',
  fileFormat: 'GeoTIFF',
  maxPixels: 1e9
});

// 2. Severity Raster Export
Export.image.toDrive({
  image: severity.toFloat(),
  description: 'Kitobo_Severity',
  folder: 'GEE_Exports',
  region: aoi.bounds(),
  scale: 10,
  crs: 'EPSG:32737',
  fileFormat: 'GeoTIFF',
  maxPixels: 1e9
});

// 3. Monthly NBR Stack Export (Named with YYYY_MM_dd prefixes for R parsing)
var nbrStack = monthly.select('NBR').toBands().toFloat().clip(aoi);
Export.image.toDrive({
  image: nbrStack,
  description: 'Kitobo_NBR_monthly',
  folder: 'GEE_Exports',
  region: aoi.bounds(),
  scale: 20,
  crs: 'EPSG:32737',
  fileFormat: 'GeoTIFF',
  maxPixels: 1e9
});
