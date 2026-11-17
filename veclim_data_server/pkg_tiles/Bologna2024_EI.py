label = "Bologna2024_EI"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..environ import DIR_DATA
from ..fun_tiles import getBologna2024Tiles

model = "Ensemble I"
mlbl = "emblI"

clbins = numpy.arange(28,38)
cllbl = [str(w) for w in clbins]

cmap = mpl.colormaps["YlOrRd_r"].resampled(len(clbins) - 1)
norm = mpl.colors.BoundaryNorm(clbins, cmap.N)

clscl = [mpl.colors.to_hex(c) for c in cmap(norm((clbins[:-1] + clbins[1:]) / 2))]

dat = numpy.load("%s/sims/Bologna2024/tile_dat_Bologna2024_annual_%s.npy" %(DIR_DATA,mlbl),mmap_mode='r')

tile_dat = {
    'label': label,
    'fun': getBologna2024Tiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}
