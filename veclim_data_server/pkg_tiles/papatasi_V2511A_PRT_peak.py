label = "papatasi_V2511A_PRT_peak"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import datetime
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getShpTiles, proj1
from ..pkg_sims import papatasi_V2511A_PRT

def calc_dat():
    dat = papatasi_V2511A_PRT.db.polys.copy()
    tms = papatasi_V2511A_PRT.peak_up_times
    dat['mean'] = [numpy.floor(tms[a][0] / 7.0) if len(tms[a]) else numpy.nan for a in tms]
    return dat.to_crs(proj1)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

clbins = numpy.arange(16,28)
cllbl = ["%d (%s)" %(a,(datetime.date(2010,1,1)+datetime.timedelta(days=int(a)*7.0)).strftime('%h')) for a in clbins]

cmap = mpl.colormaps["YlOrRd_r"].resampled(len(clbins) - 1)
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