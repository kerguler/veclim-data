import numpy

from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading papatasi_V2609C...",flush=True)

unit_scale = 2.0119863

newegg = xr_open_lazy("%s/sims/ISMED-CLIM/V2609C/sims_model_V2609C_ERA5_Europe_2010_2014_newegg_mean.nc" %(DIR_DATA))
newegg = newegg / unit_scale

newegg_ssp585 = xr_open_lazy("%s/sims/ISMED-CLIM/V2609C/sims_model_V2609C_ACCESS-ESM1-5_ssp585_Europe_2090_2094_newegg_mean.nc" %(DIR_DATA))
newegg_ssp585 = newegg_ssp585 / unit_scale

unit_scale_CERRA = 1.5210502

newegg_CERRA = xr_open_lazy("%s/sims/ISMED-CLIM/V2609C/sims_model_V2609C_Portugal_2010_2014_newegg_mean.nc" %(DIR_DATA))
newegg_CERRA = newegg_CERRA / unit_scale_CERRA
