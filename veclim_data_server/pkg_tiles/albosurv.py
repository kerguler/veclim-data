label = "albosurv"
print("Loading tiles: %s..." %label, flush=True)

import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getTiles
from ..pkg_surv import albosurv

clscl = ['#00000000', '#f15a48', '#1b3958', '#167997', '#50c0ad']
clbins = [0,1,2,3,4,5]
cllbl = ["Unknown/absent", "Global presence (2024)", "VectAbundance (2010-2022)", "AIMsurv (2020)", "VectorBase (2010-2024)"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

tran = lambda x: x

def calc_dat():
    x = albosurv.getMatrix()[:-1,:]
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
