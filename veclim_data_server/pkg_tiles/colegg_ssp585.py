ssp = "ssp585"
label = "colegg_%s" %ssp
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import get_dates, cache_npy, remove_feb29
from ..fun_tiles import getTiles
from ..pkg_sims import annualNASA

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#f47b2c', '#e85229', '#d82929', '#931b1f']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["1/16","1/8","1/4","1/2","1","2","4","8","16"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def calc_dat():
    x = annualNASA.colegg[ssp][:-1,:,:]
    return numpy.log2(numpy.nanmean(x,axis=2)/5.0) # for '2010-2019' (corrected)

def calc_date(dt):
    x = annualNASA.colegg[ssp][:-1,:,:]
    x = remove_feb29(x,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    return numpy.log2(numpy.nanmean(x,axis=2)/5.0) # for '2010-2019' (corrected)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

tile_dat = {
    'label': label,
    'fun': getTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}

def load_dates(date0,date1):
    if ((date0 == None) or (date1 == None)):
        return {}
    #
    dt = get_dates(date0, date1=date1, ts=False)
    dat_dt = cache_npy("tile_dat_%s_%s_%s.npy" %(label,dt['date0'],dt['date1']), calc_date, dt)
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