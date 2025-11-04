prd = "julsep"
label = "papatasi_%s" %prd
print("Loading tiles: %s..." %label, flush=True)

import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getShpTiles, proj1
from ..pkg_sims import papatasi2015

clscl = ["#f7fcf5", "#d5efcf", "#9ed898", "#54b567", "#1d8641", "#00441b"]
clbins = [0,1,2,3,4,5,6]
cllbl = ["0","10","100","200","300","400","500"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='neither')

def calc_dat():
    x = papatasi2015.shps[prd]
    return x.to_crs(proj1)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

tile_dat = {
    'label': label,
    'fun': getShpTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}
