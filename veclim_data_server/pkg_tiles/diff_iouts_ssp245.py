ssp = "ssp245"
label = "diff_iouts_%s" %ssp
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import get_dates, cache_npy, remove3_feb29
from ..fun_tiles import getTiles
from ..pkg_sims import annualVectorA, annualNASA

clscl = ['#50c0ad','#8dcbc1','#c6e0ee','white','#f5d9b8','#e2988a','#f15a48']
clbins = [-20.0,-10.0,-1.0,1.0,10.0,20.0]
cllbl = ['-30','-20','-10','-1','1','10','20','30']

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='both')

tran = lambda x: numpy.nanmean(x*4000,axis=2)

def calc_dat():
    x = annualVectorA.iouts[:-1,:,:].load().values
    y = annualNASA.iouts[ssp][:-1,:,:]
    return tran(y-x)

def calc_date(dt):
    x = annualVectorA.iouts[:-1,:,:].load().values
    x = remove3_feb29(x,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    y = annualNASA.iouts[ssp][:-1,:,:]
    y = remove3_feb29(y,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    return tran(y-x)

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