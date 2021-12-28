import rasterio
from rasterio.enums import Resampling
import os
from osgeo import gdal
import geopandas as gpd

#Raster dir
ras_dir=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\raster"

#Vector dir
vec_dir=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\vector"

template = gdal.Open(os.path.join(ras_dir,'Layout_Domain.tif'))
ulx, xres, xskew, uly, yskew, yres  = template.GetGeoTransform()
lrx = ulx + (template.RasterXSize * xres)
lry = uly + (template.RasterYSize * yres)
#temp_cell_size=template.transform[0]

infn = os.path.join(ras_dir,'TopoBathy_26910.tif')
outfn = os.path.join(ras_dir,'TopoBathy_26910_10m.tif')

xres=10
yres=10
resample_alg = 'average'

options = gdal.WarpOptions(options=['tr'], xRes=xres, yRes=yres,resampleAlg=resample_alg,multithread=True, outputBounds = [ulx, uly, lrx, lry])

ds = gdal.Warp(outfn, infn, options=options)
ds = None

#pts = gpd.read_file(os.path.join(vec_dir,'Template_Sampling_Points.shp'))
#pts = pts[['UTM_E', 'UTM_N', 'Value', 'geometry']]
#pts.index = range(len(pts))
#coords = [(x,y) for x, y in zip(pts.UTM_E, pts.UTM_N)]

src=rasterio.open(outfn)

# Sample the raster at every point location and store values in DataFrame
pts['Raster Value'] = [x for x in src.sample(coords)]
pts['Raster Value'] = probes.apply(lambda x: x['Raster Value'][0], axis=1)
