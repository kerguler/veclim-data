label = "Bologna2024_M0"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..environ import DIR_DATA
from ..fun_tiles import getTiles

model = 0
mlbl = 'albopictus (sPop)'

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#f47b2c', '#e85229', '#d82929', '#931b1f']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["1/16 (%s)" %(mlbl),"1/8","1/4","1/2","1","2","4","8","16"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

dat = numpy.load("%s/sims/Bologna2024/tile_dat_Bologna2024_%d_%s.npy" %(DIR_DATA,model,mlbl),mmap_mode='r')

tile_dat = {
    'label': label,
    'fun': getTiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}
