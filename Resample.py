import rasterio
from rasterio.enums import Resampling
import os
from osgeo import gdal
import geopandas as gpd

#Raster dir
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Paper"

#Vector dir
vec_dir=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\vector"

template_name="DEM_MBK_ft_NGVD29_26910.tif"

#New resolution
new_res=5

#Input name
input_nam="DEM_MBK_ft_NGVD29_26910"

no_data_ouput=999

template = gdal.Open(os.path.join(ras_dir,template_name))
ulx, xres, xskew, uly, yskew, yres  = template.GetGeoTransform()
lrx = ulx + (template.RasterXSize * xres)
lry = uly + (template.RasterYSize * yres)
#temp_cell_size=template.transform[0]

infn = os.path.join(ras_dir,input_nam+".tif")
outfn = os.path.join(ras_dir,input_nam+"_"+str(new_res)+".tif")


resample_alg = 'average'

options = gdal.WarpOptions(options=['tr'], xRes=new_res, yRes=new_res,resampleAlg=resample_alg,multithread=True,
                           outputBounds = [ulx, uly, lrx, lry],dstNodata= no_data_ouput)

ds = gdal.Warp(outfn, infn, options=options)
ds = None

