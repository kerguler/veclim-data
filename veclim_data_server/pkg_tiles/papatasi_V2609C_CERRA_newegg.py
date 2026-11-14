label = "papatasi_V2609C_CERRA_newegg"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import xarray
import matplotlib as mpl

from ..functions import cache_ncdf
from ..fun_tiles import getCERRATiles
from ..pkg_sims import papatasi_V2609C as sim_sand

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

def calc_dat():
    x2 = (
        sim_sand.newegg_CERRA
        .mean(dim='time',skipna=True)
    )
    return x2

dat = cache_ncdf("tile_dat_%s.nc" %label, calc_dat)

tile_dat = {
    'label': label,
    'fun': getCERRATiles,
    'dat': dat['newegg'],
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}