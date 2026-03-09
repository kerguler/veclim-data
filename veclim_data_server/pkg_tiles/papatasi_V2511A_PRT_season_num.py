label = "papatasi_V2511A_PRT_season_num"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import datetime
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getShpTiles, proj1
from ..pkg_sims import papatasi_V2511A_PRT

def calc_dat():
    dat = papatasi_V2511A_PRT.db.polys.copy()
    up_times = papatasi_V2511A_PRT.up_times
    dat['mean'] = [len(up_times[a]) for a in up_times]
    return dat.to_crs(proj1)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

clbins = numpy.arange(5)
cllbl = [f"{b}" for b in clbins]

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