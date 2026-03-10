label = "papatasi_V2511A_PRT_newegg"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getShpTiles, proj1
from ..pkg_sims import papatasi_V2511A_PRT

def calc_dat():
    dat = papatasi_V2511A_PRT.db.polys.copy()
    unit_scale = 1e4
    dat['mean'] = dat['Means'] / unit_scale
    return dat.to_crs(proj1)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

clbins = numpy.arange(int(dat['mean'].min()),
                      int(dat['mean'].max()))
cllbl = [f"{b}" for b in clbins[:-1]]

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