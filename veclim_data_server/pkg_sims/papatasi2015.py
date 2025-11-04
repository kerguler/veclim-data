import json
import numpy
import geopandas as gpd

from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading papatasi2015...",flush=True)

dr = "%s/sims/sand_papatasi/2015" %(DIR_DATA)
prds = ['aprdec','aprjun','julsep','octdec']
lonlat = json.load(open("%s/lonlat_sand.json" %(dr),"r"))
latitude = numpy.array(lonlat['lat'])
longitude = numpy.array(lonlat['lon'])
shps = {
    prd: gpd.read_file("%s/intersect_clc_sim_%s.shp" %(dr,prd))
    for prd in prds
}
simGERI = numpy.load("%s/mech_model_STENI_papatasi_combinedA_colnvAf_posterior_mean_GERI.npy" %(dr),mmap_mode='r')
simSTENI = numpy.load("%s/mech_model_STENI_papatasi_combinedA_colnvAf_posterior_mean_STENI.npy" %(dr),mmap_mode='r')

mask = numpy.isnan(numpy.nanmean(simSTENI[90:,:,:],axis=0))

def island(loni,lati):
    return not mask[lati,loni]