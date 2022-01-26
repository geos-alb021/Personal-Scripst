Description: This script automates the processing of satellite imagery from Sentinel-2 for the calculation and creation of rasters regarding natural color, false color and
vegetation indices (VI). The vegetation indices that are calculated are: NDVI, GNDVI, SAVI and Moisture Index.

NOTE 1: SATELLITE DATA: In the script, Sentinel-2A imagery was used. The image used was "S2A_MSIL1C_20211224T160701_N0301_R097_T16PFS_20211224T194453". It corresponds to 12/24/2021; however, you can use any other imagery, as long as it contains the AOI.
You can download data at: https://scihub.copernicus.eu/dhus/#/home

NOTE 2: Script run with python v3.8.8

NOTE 3: The following folders should be created in the same place from where you run the script:

"Sentinel2-imagery", this folder contains all the sentinel 2 bands files (where the bands files located in: "GRANULE>NameOfImage>IMG_DATA" of the sentinel2 image folder should be copied).

"clipped_raster", this folder will contain the clipped bands (temporarily).

"Shapefile_AOI", this folder contains the shapefile associated with the area of interest (AOI) for clipping the bands.

"output", this folder will contain the subfolders that will be created during the run of the script and contain the natural color, false color and vegetation index rasters generated. In every run and different Sentinel2 image, a new subfolder with the date of the imagery as its name will be created. This was made to maintain an order of the processed images.


