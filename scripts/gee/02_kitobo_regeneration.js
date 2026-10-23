/* ============================================================================
 * 02_kitobo_regeneration.js  —  Kitobo post-fire vegetation recovery
 * Earth Engine Code Editor (account ee-imnduta1). Save under users/imnduta1/Ismurn.
 *
 * Produces (a) an on-screen NBR & NDVI recovery chart over the forest, and
 * (b) a MONTHLY NBR stack exported to Drive for R-side Theil-Sen analysis
 *     (scripts/03_... Part C reads it as data/raw/gee/Kitobo_NBR_monthly.tif).
 * ========================================================================== */

// ---- 0. AOI + monitoring window --------------------------------------------
var kitobo = ee.FeatureCollection('projects/ee-imnduta1/assets/Kitobo_Forest');
var aoi    = kitobo.geometry();
Map.centerObject(aoi, 14);

var startDate = '2025-10-01';   // a season or two BEFORE the fire for a baseline
var endDate   = '2026-10-23';   // through the present

// ---- 1. Sentinel-2 SR with cloud mask + NBR/NDVI ---------------------------
function maskS2(img) {
  var scl = img.select('SCL');
  var ok = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(ok).divide(10000).copyProperties(img, ['system:time_start']);
}
function addIndices(img) {
  var nbr  = img.normalizedDifference(['B8', 'B12']).rename('NBR');
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return img.addBands([nbr, ndvi]);
}
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(aoi).filterDate(startDate, endDate)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
           .map(maskS2).map(addIndices);

// ---- 2. Monthly median composites ------------------------------------------
// Every month gets an image with the SAME bands, even when no scene survives the
// cloud filter (a fully masked placeholder), so toBands() never ends up with ragged
// bands. system:index is set to the month because toBands() uses it as the band-name
// prefix ('2025_10_01_NBR', ...), which is what R parses as the date.
// (Untested until first run in the Code Editor.)
var months = ee.List.sequence(0,
  ee.Date(endDate).difference(ee.Date(startDate), 'month').round().subtract(1));
var blank = ee.Image.constant([0, 0]).rename(['NBR', 'NDVI']).selfMask();

var monthly = ee.ImageCollection.fromImages(months.map(function(m) {
  var t0 = ee.Date(startDate).advance(m, 'month');
  var t1 = t0.advance(1, 'month');
  var col  = s2.filterDate(t0, t1).select(['NBR', 'NDVI']);
  var comp = ee.Image(ee.Algorithms.If(col.size().gt(0), col.median(), blank));
  return comp.set({'system:time_start': t0.millis(),
                   'system:index': t0.format('YYYY_MM_dd')});
}));

// ---- 3. On-screen recovery chart (mean over the forest) --------------------
var chart = ui.Chart.image.series({
  imageCollection: monthly.select(['NBR', 'NDVI']),
  region: aoi, reducer: ee.Reducer.mean(), scale: 20, xProperty: 'system:time_start'
}).setOptions({
  title: 'Kitobo Forest — monthly mean NBR & NDVI (recovery)',
  hAxis: {title: 'Date'}, vAxis: {title: 'Index'}, lineWidth: 2, pointSize: 3
});
print(chart);

// ---- 4. Export the monthly NBR stack to Drive (bands named by month) -------
// R reads the first 10 characters of each band name as the date ("%Y_%m_%d") in
// scripts/03_... Part C. toFloat() keeps all bands one type (Export needs that).
var nbrStack = monthly.select('NBR').toBands().toFloat().clip(aoi);
Export.image.toDrive({
  image: nbrStack, description: 'Kitobo_NBR_monthly', folder: 'GEE_Exports',
  region: aoi, scale: 20, crs: 'EPSG:32737', maxPixels: 1e9
});
