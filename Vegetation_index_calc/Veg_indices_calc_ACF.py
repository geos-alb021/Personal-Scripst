###-------------------Vegetation Indices Calculation For Sentinel 2 Imagery-------------------###

###Authorship: Alberto F. Coto Fonseca

#Description: This script automates the processing of satellite imagery from Sentinel-2 for the calculation and creation of
#vegetation indices (VI) rasters. Using an area of interest (AOI) shapefile all bands of the image are cropped. This new rasters are 
#saved in a new folder. Later on they are used for the calculation of the VI.

#For cropping every band with the AOI shapefile, both of them should have the same projected coordinates; therefore, a revision 
#of both crs is made. If their "crs" differ, a reprojection of the AOI is made, so it can overlap the bands.
    
#Importing libraries
import numpy as np
import rasterio as rio
import rasterio.mask as mask
from rasterio import plot
import matplotlib.pyplot as plt
import os as os
import geopandas as gpd
import fiona as fio

##Loading the bands of the satellite imagery that are going to be used in the calculation
print("Loading satellite image bands")

directory = 'C:/Users/alcof/Desktop/Scripts/Veg_index/Sentinel2_imagery'
bands = os.listdir(directory)
print('\n')
band2i=rio.open(directory + '/' + bands[1]) #For comparison of the CRS between AOI and satellite data

#Loading the shapefile for cropping
aoi = gpd.read_file('C:/Users/alcof/Desktop/Scripts/Veg_index/Shapefile_AOI/AOI.shp')
#aoi.head() #Tells the head of the geodataframe

#Checking the coordinate reference system of the AOI and imagery bands
print("AOI and bands crs:")
print("AOI crs:", aoi.crs)
print("Bands crs", band2i.crs)
print('\n')

#Converting the CRS of the AOI to the one of the imagery (if they are different)
if aoi.crs != band2i.crs:
    print("Reprojecting AOI to bands crs")
    aoi = aoi.to_crs("epsg:32616")     
print('\n')

##-------Creating a new shapefile with the reprojected AOI-------##

#Creating the new shapefile
aoi2 = "C:/Users/alcof/Desktop/Scripts/Veg_index/Shapefile_AOI/AOI2.shp"
selection = aoi[:]
selection.to_file(aoi2)

#Opening the AOI shapefile with the new projection coordinates
with fio.open("C:/Users/alcof/Desktop/Scripts/Veg_index/Shapefile_AOI/AOI2.shp", "r") as shapefile:
    shape = [feature["geometry"] for feature in shapefile]

#Cropping the imagery package to an area of interest (AOI)
print("Cropping the bands...")
for band in bands:
    
    rasterPath = os.path.join(directory,band)
    rasterBand = rio.open(rasterPath)
    outImage, outTransform = rio.mask.mask(rasterBand, shape, crop=True)
    outMeta = rasterBand.meta
    
    outMeta.update({"driver": 'JP2OpenJPEG',
                "height": outImage.shape[1],
                "width": outImage.shape[2],
                "transform": outTransform})
    
    outPath = os.path.join('C:/Users/alcof/Desktop/Scripts/Veg_index/clipped_raster',band)
    outRaster = rio.open(outPath, "w", **outMeta) 
    outRaster.write(outImage)
    outRaster.close()
print("Bands cropped")
print('\n')

clipped_directory = 'C:/Users/alcof/Desktop/Scripts/Veg_index/clipped_raster'
clipped_bands = os.listdir(clipped_directory)

band2 = rio.open(clipped_directory + '/' + clipped_bands[1]) #Blue
band3 = rio.open(clipped_directory + '/' + clipped_bands[2]) #Green
band4 = rio.open(clipped_directory + '/' + clipped_bands[3]) #Red
band8 = rio.open(clipped_directory + '/' + clipped_bands[7]) #NIR
band11 = rio.open(clipped_directory + '/' + clipped_bands[10]) #SWIR
band8A = rio.open(clipped_directory + '/' + clipped_bands[12]) #VNIR
    
#Obtaining information of the cropped imagery
print("Raster information")
height=band2.height
width=band2.width
print("Bands height:{}".format(height))
print("Bands widht:{}".format(width))
print('\n')

#Converting the bands from integer to float data ( float arrays)
print("Converting bands to float data")
red= band4.read(1).astype('float64')
blue = band2.read(1).astype('float64')
green = band3.read(1).astype('float64')
nir= band8.read(1).astype('float64')
swir = band11.read(1).astype('float64')
vnir = band8A.read(1).astype('float64')
print('\n')

#Allowing division by zero
np.seterr(divide='ignore', invalid='ignore')


###--------------Calculating Spectral Indices-----------------###
print("Calculating spectral indices")

ndvi = np.where(nir+red==0, 0, (nir-red)/(nir+red)) #Normalized Difference Vegetation Index (NDVI)

gndvi = np.where(nir+green==0, 0, (nir-green)/(nir+green)) #Green Normalized Difference Vegetation Index (GNDVI)

L = 0.1 #L is a constant to reduce bias due to ground background effects. Value of 0.1 is used according to literature.
#Reference: M. Tasumi, “Progress in operational estimation of regional evapotranspiration using satellite imagery”,University of Idaho, 2003.

savi = np.where(L+(nir+red)==0, -9999,(((1+1)*(nir-red))/(L+(nir+red)))) #Adjusted soil vegetation index (SAVI)

mi = np.where (vnir+swir==0, 0, (vnir-swir)/(vnir+swir)) #Moisture Index (MI)
print('\n')

#Exploring the indices data (this helps to look for outliers)
print("Exploring indices data")
print("NDVI min:",ndvi.min())
print("NDVI max:",ndvi.max())
print('\n')
print("GNDVI min:",gndvi.min())
print("GNDVI max:",gndvi.max())
print('\n')
print("SAVI min:",savi.min())
print("SAVI max:",savi.max())
print('\n')
print("Moisture index min:", mi.min())
print("Moisture index max:", mi.max())
print('\n')

###---------Creating the combination of bands rasters and vegetation index rasters---------##

print("Creating rasters of the band combinationa...")
#Creating a natural color combination raster
naturalColor = rio.open('C:/Users/alcof/Desktop/Scripts/Veg_index/output/sent2natcolor.tiff', 'w', driver='Gtiff',
                     width=width, height=height,
                     count=3, #number of bands
                     crs=band2.crs,
                     transform=band2.transform,
                     dtype=band2.dtypes[0]
                     )
#True color combination: red, green and blue.
naturalColor.write(band4.read(1),1) #Red
naturalColor.write(band3.read(1),2) #Green
naturalColor.write(band2.read(1),3) #Blue
naturalColor.close()

#Creating a false color combination raster
falseColor = rio.open('C:/Users/alcof/Desktop/Scripts/Veg_index/output/sent2falsecolor.tiff', 'w', driver='Gtiff',
                     width=width, height=height,
                     count=3, #number of bands
                     crs=band2.crs,
                     transform=band2.transform,
                     dtype=band2.dtypes[0]
                     )
#False color combination: NIR, Red and Green.
falseColor.write(band8.read(1),1) #NIR
falseColor.write(band4.read(1),2) #Red
falseColor.write(band3.read(1),3) #Green
falseColor.close()
print('\n')

print("Creating rasters of the vegetation indices...")
#Creating the NDVI raster
NDVI = rio.open('C:/Users/alcof/Desktop/Scripts/Veg_index/output/sent2NDVI.tiff', 'w', driver='Gtiff',
                     width=width, height=height,
                     count=1, #number of bands
                     crs=band2.crs,
                     transform=band2.transform,
                     dtype='float64'
                     )
NDVI.write(ndvi,1)
NDVI.close()

#Creating the GNDVI raster
GNDVI = rio.open('C:/Users/alcof/Desktop/Scripts/Veg_index/output/sent2GNDVI.tiff', 'w', driver='Gtiff',
                     width=width, height=height,
                     count=1, #number of bands
                     crs=band2.crs,
                     transform=band2.transform,
                     dtype='float64'
                     )
GNDVI.write(gndvi,1)
GNDVI.close()

#Creating the SAVI raster
SAVI = rio.open('C:/Users/alcof/Desktop/Scripts/Veg_index/output/sent2SAVI.tiff', 'w', driver='Gtiff',
                     width=width, height=height,
                     count=1, #number of bands
                     crs=band2.crs,
                     transform=band2.transform,
                     dtype='float64'
                     )
SAVI.write(savi,1)
SAVI.close()

#Creating the MI raster
MI = rio.open('C:/Users/alcof/Desktop/Scripts/Veg_index/output/sent2MI.tiff', 'w', driver='Gtiff',
                     width=width, height=height,
                     count=1, #number of bands
                     crs=band2.crs,
                     transform=band2.transform,
                     dtype='float64'
                     )
MI.write(savi,1)
MI.close()
print('\n')

###--------------Plotting the results--------------###

#Opening the rasters
#ndvi_image = rio.open('C:/Users/alcof/Documents/DOCTORADO/PROGRAMMING/Veg_index/output/sent2NDVI.tiff', count=1)
#gndvi_image = rio.open('C:/Users/alcof/Documents/DOCTORADO/PROGRAMMING/Veg_index/output/sent2GNDVI.tiff', count=1)
#savi_image = rio.open('C:/Users/alcof/Documents/DOCTORADO/PROGRAMMING/Veg_index/output/sent2SAVI.tiff', count=1)


#Normalizing and stacking the bands of the true color and false color rasters for the correct color visualisation

def normalize(array): #Utilizing a normalizing funtion
    array_min = array.min()
    array_max = array.max()
    return (array-array_min)/(array_max-array_min)

red2 = normalize(red)   
blue2 = normalize(blue)
green2 = normalize(green)
nir2 = normalize(nir)

#Stacking the bands for visualization
nc_norm = np.dstack((red2, green2, blue2))
fc_norm = np.dstack((nir2, red2, green2))

print("Plotting the combination bands and the vegetation indices rasters")
fig, (ax1,ax2) = plt.subplots(1,2,figsize= (21,7))
fig.suptitle('Combination bands')
ax1.imshow(nc_norm)
ax2.imshow(fc_norm)

ax1.set_title("Natural Color")
ax2.set_title("False Color")
fig.tight_layout()  

fig, ((ax3,ax4),(ax5,ax6)) = plt.subplots( nrows=2, ncols=2, figsize= (21,7))
fig.suptitle('Combination bands and Vegetation Spectral Indices')
plot.show(ndvi, cmap='RdYlGn', ax=ax3,)
plot.show(savi, cmap='Spectral', ax=ax4)
plot.show(gndvi, cmap='RdBu', ax=ax5,)
plot.show(mi, cmap='bwr', ax=ax6,)

ax3.set_title("NDVI")
ax4.set_title("SAVI")
ax5.set_title("Green NDVI")
ax6.set_title("Moisture Index")

fig.tight_layout()               







