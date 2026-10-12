// ============================================================
// Functional Biomes and Stable Forest Mask
// ============================================================
//
// Creates:
//   1. A stable forest mask for 2001–2020
//   2. Tree-cover percentage at 500 m for 2000 and 2020
//
// Data sources:
//   - MODIS MCD12Q1.061
//   - GLAD GLCLU2020
//   - Functional biome polygons
//
// NOTE:
// Replace YOUR_PROJECT_ID with an Earth Engine Cloud Project
// to which you have write access before running this export.
//
// ============================================================

var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);           \
\
//Load functional biome layer\
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ // \
// Replace with the path to your own Earth Engine asset containing the sample points.
var polygons = ee.FeatureCollection("projects/YOUR_PROJECT_ID/assets/functional_biomes_1_12")\
  .select('grps90'); //Polygons of functional biomes\

\
/////Visualisaton on the map/////\
// Define a palette for the 12 distinct functional biomes\
var palette = [\
  "#cba10b", // "0"\
  "#010eff", // "1"\
  "#ff3c01", // "9"\
  "#d100e4", // "12"\
  "#308139", // "13"\
  "#e66f89", // "15"\
  "#ac1620", // "17"\
  "#ff0000", // "18"\
  "#fff201", // "19"\
  "#ffd501", // "20"\
  "#46a051", // "23"\
  "#57db6f"  // "24"\
];\
\
\
// Create an empty image into which to paint the features, cast to byte.\
var empty = ee.Image().byte();\
// Paint the interior of the polygons with different colors.\
var fills = empty.paint(\{\
  featureCollection: polygons,\
  color: 'grps90',\
\});\
Map.addLayer(fills, \{palette: palette, min:0, max: 24, opacity: 0.7\}, 'Functional Biomes'), false;\
\
//Turn polygons features into image\
var Funct_biomes = polygons.reduceToImage(\{\
  properties: ['grps90'],\
  reducer: ee.Reducer.first()\
\}).rename('grps90')\
\
print(Funct_biomes);\
\
// CREATE MASK OF STABLE FORESTLAND AND TREE COVER\
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ // \
// 1. Create image of forest land cover for each year (binary)\
// Load the MODIS Land Cover Type dataset (MCD12Q1.061)g\
\
// Create a collection of land cover images for all years\
var modisLandCover = ee.ImageCollection("MODIS/061/MCD12Q1")\
                       .select('LC_Type1');\
\
\
// 2. Create image of stable forestland (binary)\
// Define pixels with values 1-5 and 8-9 as forestland, everything else is 0\
var years = ee.List.sequence(2001, 2020);\
\
function getForestMask(year) \{\
  var image = modisLandCover\
    .filter(ee.Filter.calendarRange(year, year, 'year'))\
    .first();\
  \
  // Forest classes in MODIS IGBP scheme: 1-5 and 8-9 (forest and savanna)\
  var forestMask = image.remap([1,2,3,4,5,8,9], [1,1,1,1,1,1,1], 0);\
  return forestMask\
          .rename('forestMask')\
          .set('year', year);\
\}\
\
var forestMasks = ee.ImageCollection.fromImages(years.map(getForestMask));\
\
print(forestMasks)\
\
// Map first year\
var firstForestImage = ee.Image(forestMasks.first());\
Map.addLayer(firstForestImage.updateMask(firstForestImage), \
             \{min: 0, max: 1, palette: ['white', 'green']\}, \
             'Forest Mask - First Year'), false;\
\
// Map last year\
// Get the last image by sorting descending on 'year'\
var lastForestImage = forestMasks\
  .sort('year', false)   // false = descending\
  .first();\
/*\
Map.addLayer(\
  lastForestImage.updateMask(lastForestImage), \
  \{min: 0, max: 1, palette: ['white', 'darkgreen']\}, \
  'Forest Mask - Last Year'\
);*/\
             \
// Get the first year (2001) and use it as a reference\
var referenceImage = forestMasks.first();\
\
var stableMask = forestMasks.toBands().reduce(\{\
  reducer: ee.Reducer.product()\
\}).rename('stableMask')\
\
// Unstable areas are where the land cover has changed over time\
var unstableMask = stableMask.not().rename('unstableMask');\
\
             \
       \
// Combine stable and unstable areas into one image\
var stabilityImage = stableMask.addBands(unstableMask);\
\
// Add stable areas in blue\
Map.addLayer(\
  stableMask,\
  \{min: 0, max: 1, palette: ['white', 'blue']\},\
  'Stable Forest', false\
);\
\
// Add unstable areas in red\
Map.addLayer(\
  unstableMask,\
  \{min: 0, max: 1, palette: ['white', 'red']\},\
  'Unstable Forest', false\
);\
\
//Export mask\
Export.image.toAsset(\{\
  image: stableMask,\
  description: 'stableMask',\
  assetId: 'projects/YOUR_PROJECT_ID/assets/stableMask_LC_type',\
  scale: 500,\
  maxPixels: 1e11,\
  region: unboundedGeo\
\});\
\
\
// 4. Create stable tree cover image for 2000 and 2020 (binary)\
// Load tree cover data\
var m00 = ee.Image('projects/glad/GLCLU2020/LCLUC_2000');\
var m20 = ee.Image('projects/glad/GLCLU2020/LCLUC_2020');\
\
// Define threshold values for the range of stable forest cover (from land cover GLAD data) to be included\
var lowerThreshold1 = 25; // Lower limit for range (stable tree cover 3m trees)\
var upperThreshold1 = 48; // Upper limit for range (stable tree cover >25m trees)\
\
// Create a binary image based on the threshold range\
var binaryImage1_2000 = m00.gte(lowerThreshold1).and(m00.lte(upperThreshold1));\
var binaryImage1_2020 = m20.gte(lowerThreshold1).and(m20.lte(upperThreshold1));\
\
\
print(binaryImage1_2000);\
\
// Map\
Map.addLayer(\
  binaryImage1_2000.updateMask(binaryImage1_2000).unmask(0),\
  \{min: 0, max: 1, palette: ['white', 'darkgreen']\}, \
  'Tree cover 2000', false\
);\
\
Map.addLayer(\
  binaryImage1_2020.updateMask(binaryImage1_2020).unmask(0), \
  \{min: 0, max: 1, palette: ['white', 'darkgreen']\}, \
  'Tree cover 2020', false\
);\
\
// Calculate the fraction of stable tree cover (30m) within a 500m pixel\
// Aggregate using reduceResolution\
\
// For 2000\
var forestPerc_2000 = binaryImage1_2000\
  .reduceResolution(\{reducer: ee.Reducer.mean(), maxPixels: 323\})\
  .reproject(\{crs: 'EPSG:4326', scale: 500\})\
  .multiply(100).rename('TreeCoverPct'); // Convert from fraction to percent\
\
print(forestPerc_2000)\
\
Export.image.toAsset(\{\
  image: forestPerc_2000,\
  description: 'forestPerc_2000',\
  assetId: 'projects/YOUR_PROJECT_ID/assets/forestPerc_2000',\
  crs: 'EPSG:4326',\
  scale: 500,\
  maxPixels: 1e11,\
  region: unboundedGeo\
\});\
\
// For 2020\
var forestPerc_2020 = binaryImage1_2020\
  .reduceResolution(\{reducer: ee.Reducer.mean(), maxPixels: 323\})\
  .reproject(\{crs: 'EPSG:4326', scale: 500\})\
  .multiply(100).rename('TreeCoverPct'); // Convert from fraction to percent\
\
\
Export.image.toAsset(\{\
  image: forestPerc_2020,\
  description: 'forestPerc_2020',\
  assetId: 'projects/YOUR_PROJECT_ID/assets/forestPerc_2020',\
  crs: 'EPSG:4326',\
  scale: 500,\
  maxPixels: 1e11,\
  region: unboundedGeo\
\});\
\
\
Map.addLayer(forestPerc_2000, \{min:0, max:100\}, "forestPerc_2000", false)\
Map.addLayer(forestPerc_2020, \{min:0, max:100\}, "forestPerc_2020", false)\
\