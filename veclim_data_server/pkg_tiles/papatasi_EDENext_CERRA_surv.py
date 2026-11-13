from ..environ import DIR_DATA

label = "papatasi_EDENext_CERRA_surv"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

import geopandas as gpd

from ..functions import cache_npy
from ..fun_tiles import getShpTiles, proj1
from ..pkg_surv import sandfly_EDENext as surv

def calc_dat():
    srv = surv.surv
    dat = surv.getMatrixCERRA(srv)
    dat['mean'] = 1.0
    return dat.to_crs(proj1)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

clbins = numpy.array([0,1,2])
cllbl = ["Unknown/Absent", "Reported (EDENext)"]

cmap = mpl.colormaps["YlOrRd"].resampled(len(clbins) - 1)
norm = mpl.colors.BoundaryNorm(clbins, cmap.N)

clscl = [mpl.colors.to_hex(c) for c in cmap(norm((clbins[:-1] + clbins[1:]) / 2))]

tile_dat = {
    'label': label,
    'fun': getShpTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}