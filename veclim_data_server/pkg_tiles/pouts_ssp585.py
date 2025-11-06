ssp = "ssp585"
label = "pouts_%s" %ssp
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import get_dates, cache_npy, remove3_feb29
from ..fun_tiles import getTiles
from ..pkg_sims import annualNASA

clscl = ['#00000000', '#fbe590', '#f7a034', '#e85229', '#931b1f']
clbins = [0,0.01,0.05,0.1,0.2,0.5]
cllbl = ["0","1","5","10","20","50"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

tran = lambda x: numpy.nanmean(x,axis=2)

def calc_dat():
    x = annualNASA.pouts[ssp][:-1,:,:]
    return tran(x)

def calc_date(dt):
    x = annualNASA.pouts[ssp][:-1,:,:]
    x = remove3_feb29(x,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    return tran(x)

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
        'label': "%s_%s_%s" %(label,dt['date0'],dt['date1']),
        'cllbl': cllbl,
        'clscl': clscl
    }