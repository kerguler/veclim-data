label = "papatasi_V2609C"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import xarray
import matplotlib as mpl

from ..functions import cache_ncdf
from ..fun_tiles import getTiles
from ..pkg_sims import papatasi_V2609C

clscl = ['#00000000', 
         '#fbe590', 
         '#f7a034', 
         '#e85229', 
         '#931b1f']
clbins = [0,
          1,
          2,
          3,
          4,
          5]
cllbl = ["0",
         "1",
         "2",
         "3",
         "4",
         "5"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def calc_dat():
    x2 = (
        papatasi_V2609C.newegg
        .isel(lat=slice(None, -1))
        .mean(dim='time',skipna=True)
    )
    return x2

dat = cache_ncdf("tile_dat_%s.nc" %label, calc_dat)

tile_dat = {
    'label': label,
    'fun': getTiles,
    'dat': dat['newegg'],
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}