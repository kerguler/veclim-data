label = "larva"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import get_dates, cache_npy, remove3_feb29
from ..fun_tiles import getTiles
from ..pkg_sims import annualVectorA

clscl = numpy.array([
         '#081d58',
         '#253494',
         '#225ea8',
         '#1d91c0',
         '#41b6c4',
         '#7fcdbb',
         '#fcc65a', 
         '#f7a034', 
         '#f47b2c', 
         '#e85229', 
         '#d82929',
         '#b42125',
         '#00000000',
         ])[::-1].tolist()
clbins = numpy.cumsum([0,1,30,28,31,30,31,30,31,31,30,31,30,31])
cllbl = ["NA","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def trann2(mat):
    tmp = mat[:,:,:]>1.0
    annmat = numpy.sum(numpy.cumsum(tmp[:,:,1:]>tmp[:,:,:-1],axis=2)==0,axis=2)+1
    annmat[annmat==365] = 0
    return annmat

def calc_dat():
    x = annualVectorA.coln2[:-1,:,:]
    return trann2(x)

def calc_date(dt):
    x = annualVectorA.coln2[:-1,:,:]
    x = remove3_feb29(x,
                     dt['days'],
                     dt['isFeb29'],
                     tolist=False)
    return trann2(x)

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