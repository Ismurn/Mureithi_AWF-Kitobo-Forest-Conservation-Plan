/* ============================================================================
 * 01_kitobo_dNBR.js  —  Kitobo Forest 2026 burn severity (Sentinel-2 dNBR)
 * Earth Engine Code Editor (account ee-imnduta1). Save under users/imnduta1/Ismurn.
 *
 * WORKFLOW
 *   1. In R, run scripts/01_boundary.R -> outputs/shapefiles/Kitobo_Forest.shp
 *   2. Upload that shapefile as an EE table asset:
 *        projects/ee-imnduta1/assets/Kitobo_Forest
 *   3. Set the PRE/POST windows from the FIRMS fire date (scripts/03_... Part A).
 *   4. Run, then download both Drive exports into the repo's data/raw/gee/ folder.
 *
 * METHOD  Key & Benson dNBR. NBR = (NIR - SWIR2)/(NIR + SWIR2) = (B8 - B12)/(B8 + B12).
 *         dNBR = NBR_pre - NBR_post (scaled x1000). Positive = vegetation loss.
 * ========================================================================== */

// ---- 0. Area of interest ---------------------------------------------------
var kitobo = ee.FeatureCollection('projects/ee-imnduta1/assets/Kitobo_Forest');
var aoi    = kitobo.geometry().buffer(2000);   // small buffer to see the burn edge
Map.centerObject(aoi, 14);

// ---- 1. Fire-bracketing windows (EDIT from the FIRMS date) -----------------
var preStart  = '2026-02-01', preEnd  = '2026-02-28';   // clear window BEFORE the fire
var postStart = '2026-04-01', postEnd = '2026-04-30';   // clear window AFTER the fire

// ---- 2. Sentinel-2 SR: cloud/shadow mask via Scene Classification (SCL) ----
function maskS2(img) {
  var scl = img.select('SCL');
  // drop: 3 shadow, 8/9 cloud med/high, 10 cirrus, 11 snow
  var ok = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(ok).divide(10000)
            .copyProperties(img, ['system:time_start']);
}
function addNBR(img) {
  return img.addBands(img.normalizedDifference(['B8', 'B12']).rename('NBR'));
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(aoi)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
           .map(maskS2).map(addNBR);

var preNBR  = s2.filterDate(preStart,  preEnd ).select('NBR').median();
var postNBR = s2.filterDate(postStart, postEnd).select('NBR').median();

// ---- 3. dNBR + Key & Benson severity classes -------------------------------
var dNBR = preNBR.subtract(postNBR).multiply(1000).rename('dNBR').clip(aoi);

var severity = ee.Image(0)
  .where(dNBR.gt(100).and(dNBR.lte(270)), 1)   // low
  .where(dNBR.gt(270).and(dNBR.lte(440)), 2)   // moderate-low
  .where(dNBR.gt(440).and(dNBR.lte(660)), 3)   // moderate-high
  .where(dNBR.gt(660),                    4)   // high
  .rename('severity').clip(aoi).toByte();

// ---- 4. Display ------------------------------------------------------------
var sevPalette = ['d9d9d9', 'ffffb2', 'fecc5c', 'fd8d3c', 'e31a1c'];
Map.addLayer(dNBR, {min: -200, max: 800, palette: ['006400','ffffff','ff0000']}, 'dNBR');
Map.addLayer(severity, {min: 0, max: 4, palette: sevPalette}, 'Severity');
Map.addLayer(kitobo.style({color: '000000', fillColor: '00000000', width: 2}), {}, 'Kitobo');

// ---- 5. Export both rasters to Drive (download into repo data/raw/gee/) ----
Export.image.toDrive({
  image: dNBR, description: 'Kitobo_dNBR', folder: 'GEE_Exports',
  region: aoi, scale: 10, crs: 'EPSG:32737', maxPixels: 1e9
});
Export.image.toDrive({
  image: severity, description: 'Kitobo_Severity', folder: 'GEE_Exports',
  region: aoi, scale: 10, crs: 'EPSG:32737', maxPixels: 1e9
});

// Optional: area burned (ha) by severity class, printed to the Console
var areaImg = ee.Image.pixelArea().divide(1e4).addBands(severity);
var areas = areaImg.reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'severity'}),
  geometry: kitobo.geometry(), scale: 10, maxPixels: 1e9
});
print('Burned area (ha) by severity class within Kitobo:', areas);
