label = "pop"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl
from matplotlib import pyplot as plt

from ..functions import cache_npy
from ..fun_tiles import getTiles
from ..pkg_sims import popdens

clscl = []
clbins = []
cllbl = []

cmap = plt.cm.viridis
norm = mpl.colors.Normalize(vmin=-3,vmax=3)

tran = lambda x: numpy.log10(x)

def calc_dat():
    x = popdens.pop[:-1,:]
    return tran(x)

dat = cache_npy("tile_dat_%s.npy" %label, calc_dat)

tile_dat = {
    'label': label,
    'fun': getTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}
