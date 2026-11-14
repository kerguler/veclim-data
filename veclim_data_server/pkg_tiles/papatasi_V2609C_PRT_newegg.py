label = "papatasi_V2609C_PRT_newegg"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getShpTiles, proj1
from ..pkg_sims import papatasi_V2609C_PRT

def calc_dat():
    dat = papatasi_V2609C_PRT.db.polys.copy()
    unit_scale = 1.5210502
    dat['mean'] = dat['Means'] / unit_scale
    return dat.to_crs(proj1)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

clscl = ['#00000000', 
         '#fbe590', 
         '#f7a034', 
         '#e85229', 
         '#931b1f']
clbins = [0,
          0.5,
          1,
          1.5,
          2,
          2.5]
cllbl = ["0",
         "0.5",
         "1",
         "1.5",
         "2",
         "2.5"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

tile_dat = {
    'label': label,
    'fun': getShpTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}