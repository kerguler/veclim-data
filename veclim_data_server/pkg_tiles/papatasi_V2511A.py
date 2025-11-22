label = "papatasi_V2511A"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import xarray
import matplotlib as mpl

from ..functions import cache_ncdf
from ..fun_tiles import getCERRATiles
from ..pkg_sims import papatasi_V2511A

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#f47b2c', '#e85229', '#d82929', '#931b1f']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["1/16","1/8","1/4","1/2","1","2","4","8","16"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def calc_dat():
    x = papatasi_V2511A.sand.mean(dim='time',skipna=True)
    x2 = xarray.apply_ufunc(
        lambda z: numpy.log2(z / 10000.0),
        x,
        dask="parallelized",
        vectorize=True,
    )
    return x2

dat = cache_ncdf("tile_dat_%s.nc" %label, calc_dat)

tile_dat = {
    'label': label,
    'fun': getCERRATiles,
    'dat': dat['female_mn'],
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}