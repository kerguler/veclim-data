label = "Bologna2024_vabun"
print("Loading tiles: %s..." %label, flush=True)

import numpy
import matplotlib as mpl

from ..functions import cache_npy
from ..fun_tiles import getBologna2024Tiles
from ..pkg_surv import Bologna2024

clscl = ['#00000000', '#fbe590', '#f7a034', '#e85229', '#931b1f']
clbins = [0,1,2,3,4,5]
cllbl = ["Unknown/absent", "VectAbundance (2019-2022) - 1 year", "VectAbundance (2019-2022) - 2 years", "VectAbundance (2019-2022) - 3 years", "VectAbundance (2019-2022) - 4 years"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')

def calc_dat():
    mat = Bologna2024.surv
    years = numpy.array([s.split('-')[0] for s in mat.yrwk.values])
    unique_years = numpy.unique(years)
    #
    counts = numpy.zeros((mat.sizes['lat'], mat.sizes['lon']), dtype=int)
    for y in unique_years:
        mask = years == y
        any_non_nan = mat.eggs.isel(yrwk=mask).notnull().any(dim='yrwk')
        counts += any_non_nan.values.astype(int)
    #
    return counts

dat = cache_npy("tile_dat_Bologna2024.npy", calc_dat)

tile_dat = {
    'label': label,
    'fun': getBologna2024Tiles,
    'dat': dat,
    'cmap': cmap,
    'norm': norm,
    'cllbl': cllbl,
    'clscl': clscl
}
