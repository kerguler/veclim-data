label = "chikv_iouts_fcast"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import pandas
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getTiles
from ..pkg_sims import forecastECMWF

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#e85229', '#931b1f']
clbins = [0,1,10,50,100,500,1000]
cllbl = ["0","1","10","50","100","500","1000"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

tran = lambda x: numpy.nanmean(x*4000,axis=2)

def calc_date(dt):
    x = forecastECMWF.iouts[:-1,:,dt].load().values
    return tran(x)

tile_dat = {}

def load_dates(date0,date1):
    if ((date0 == None) or (date1 == None)):
        return {}
    #
    dt0 = pandas.to_datetime(date0)
    dt1 = pandas.to_datetime(date1)
    #
    dt = (forecastECMWF.idates >= dt0) & (forecastECMWF.idates <= dt1)
    if numpy.abs(numpy.sum(dt) - (dt1-dt0).days) > 7:
        return {'error': "Forecast dates do not match the request!"}
    #
    dt0lab = dt0.strftime("%Y-%m-%d")
    dt1lab = dt1.strftime("%Y-%m-%d")
    #
    dat_dt = cache_npy("tile_dat_%s_%s_%s.npy" %(label,dt0lab,dt1lab), calc_date, dt)
    #
    return {
        'fun': getTiles,
        'dat': dat_dt,
        'cmap': cmap,
        'norm': norm,
        'label': "%s_%s_%s" %(label,dt0lab,dt1lab),
        'cllbl': cllbl,
        'clscl': clscl
    }