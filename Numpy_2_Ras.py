from osgeo import gdal
import numpy as np
import os

#List of numpy arrays we have to process
npdir1=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223\np_arrays"
names1=[r"\Draft_1_1_ft\Draft_1_1_ft_2056",r'\Draft_2_6_ft\Draft_2_6_ft_2056']
Paths=[npdir1+"\\"+name  for name in names1]


#Directory of template
ras_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\raster'



#Raster template. It needs to have the same dimentions as the numpy matrix
ras_temp_nam='Layout_Domain.tif'
#Let's open the template raster in gdal
data0=gdal.Open(os.path.join(ras_dir,ras_temp_nam))
#Transformation
geotrans=data0.GetGeoTransform()
#Projection
proj=data0.GetProjection()
#Raster dimension
x_pixels,y_pixels=data0.RasterXSize,data0.RasterYSize
#Tiff driver
driver = gdal.GetDriverByName('GTiff')

#Now, let's loop through the numpy arrays and exporter to rasters with the same dimensions as the template
for path in Paths:
    #Names of arrays
    nam=path[path.rfind("\\")+1:]
    #path where we'll export the new raster
    dst_filename = os.path.join(ras_dir,nam+".tif")
    # Let's create the empty raster
    dataset = None
    dataset = driver.Create(dst_filename,x_pixels, y_pixels, 1,gdal.GDT_Float32)
    #Let's load the numpy array we want to export
    array=np.load(path+".npy")
    # Let's write the content of the numpy array in the new raster
    dataset.GetRasterBand(1).WriteArray(array)
    #Let's set the transform of the new raster
    dataset.SetGeoTransform(geotrans)
    # Let's set the projection of the new raster
    dataset.SetProjection(proj)
    #Writes new raster to disk
    dataset.FlushCache()
    #Closes the dataset
    dataset=None

