label = "presence"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getTiles
from ..pkg_surv import albosurv

clscl = ['#00000000', '#931b1f']
clbins = [0,0.5,1]
cllbl = ["Unknown/absent","Reported/established"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

tran = lambda x: numpy.array(x)

def calc_dat():
    x = albosurv.presence.matrix[:-1,:]
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
