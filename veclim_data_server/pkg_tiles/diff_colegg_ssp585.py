ssp = "ssp585"
label = "diff_colegg_%s" %ssp
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import get_dates, cache_npy, remove3_feb29
from ..fun_tiles import getTiles
from ..pkg_sims import annualVectorA, annualNASA

clscl = ['#50c0ad','#8dcbc1','#c6e0ee','white','#f5d9b8','#e2988a','#f15a48']
clbins = [-3,-2,-1,1,2,3]
cllbl = ["1/16","1/8","1/4","1/2","2","4","8","16"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='both')

def calc_dat():
    x = annualNASA.colegg[ssp][:-1,:,:]
    tmp = numpy.log2(numpy.nanmean(x,axis=2)/5.0)
    tmp[tmp<-4] = -4
    tmp[tmp>4] = 4
    y = annualVectorA.colegg[:-1,:,:].load().values
    tmpf = numpy.log2(numpy.nanmean(y,axis=2)/5.0)
    tmpf[tmpf<-4] = -4
    tmpf[tmpf>4] = 4
    tmp = tmp-tmpf
    return tmp

def calc_date(dt):
    x = annualNASA.colegg[ssp][:-1,:,:]
    x = remove3_feb29(x,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    tmp = numpy.log2(numpy.nanmean(x,axis=2)/5.0)
    tmp[tmp<-4] = -4
    tmp[tmp>4] = 4
    y = annualVectorA.colegg[:-1,:,:].load().values
    y = remove3_feb29(y,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    tmpf = numpy.log2(numpy.nanmean(y,axis=2)/5.0)
    tmpf[tmpf<-4] = -4
    tmpf[tmpf>4] = 4
    tmp = tmp-tmpf
    return tmp

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