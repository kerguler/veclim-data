label = "colegg_fcast"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import pandas
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getTiles
from ..pkg_sims import forecastECMWF

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#f47b2c', '#e85229', '#d82929', '#931b1f']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["1/16","1/8","1/4","1/2","1","2","4","8","16"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def calc_date(dt):
    x = forecastECMWF.colegg[:-1,:,dt].load().values
    return numpy.log2(numpy.nanmean(x,axis=2)/5.0) # for '2010-2019' (corrected)

tile_dat = {}

def load_dates(date0,date1):
    if ((date0 == None) or (date1 == None)):
        return {}
    #
    dt0 = pandas.to_datetime(date0)
    dt1 = pandas.to_datetime(date1)
    #
    dt = (forecastECMWF.dates >= dt0) & (forecastECMWF.dates <= dt1)
    if numpy.abs(numpy.sum(dt) - (dt1-dt0).days) > 7:
        return {'error': "Forecast dates do not match the request!"}
    #
    print(date0,date1,dt0,dt1)
    #
    dat_dt = cache_npy("tile_dat_%s_%s_%s.npy" %(label,date0,date1), calc_date, dt)
    #
    return {
        'fun': getTiles,
        'dat': dat_dt,
        'cmap': cmap,
        'norm': norm,
        'label': '',
        'cllbl': cllbl,
        'clscl': clscl
    }