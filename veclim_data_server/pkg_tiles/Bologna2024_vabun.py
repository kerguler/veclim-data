label = "Bologna2024_vabun"
print("Loading tiles: %s..." %label, flush=True)

import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getTiles
from ..pkg_surv import Bologna2024

clscl = ['#00000000', '#f15a48']
clbins = [0,1]
cllbl = ["Unknown/absent", "VectAbundance (2010-2022)"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def calc_dat():
    x = Bologna2024.surv.eggs.load().values > 0.0
    return x

dat = cache_npy("tile_dat_Bologna2024.npy", calc_dat)

tile_dat = {
    'label': label,
    'fun': getTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}
