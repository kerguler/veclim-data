label = "diff_pouts_1980"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import get_dates, cache_npy, remove_feb29
from ..fun_tiles import getTiles
from ..pkg_sims import annualVectorA, annualVectorA_1980

clscl = ['#50c0ad','#8dcbc1','#c6e0ee','white','#f5d9b8','#e2988a','#f15a48']
clbins = [-0.15,-0.1,-0.05,0.05,0.1,0.15]
cllbl = ['-20%','-15%','-10%','-5%','5%','10%','15%','20%']

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='both')

tran = lambda x: numpy.nanmean(x,axis=2)

def calc_dat():
    x = annualVectorA.pouts[:-1,:,:].load().values
    y = annualVectorA_1980.pouts[:-1,:,:].load().values
    return tran(x-y)

def calc_date(dt):
    x = annualVectorA.pouts[:-1,:,:].load().values
    x = remove_feb29(x,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    y = annualVectorA_1980.pouts[:-1,:,:].load().values
    y = remove_feb29(y,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    return tran(x-y)

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