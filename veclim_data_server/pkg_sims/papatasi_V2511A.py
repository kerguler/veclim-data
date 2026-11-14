from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading papatasi_V2511A...",flush=True)

female_md = xr_open_lazy("%s/sims/ISMED-CLIM/V2511A/sims_model_V2511A_Portugal_female_md.nc" %(DIR_DATA))
female_hi = xr_open_lazy("%s/sims/ISMED-CLIM/V2511A/sims_model_V2511A_Portugal_female_hi.nc" %(DIR_DATA))
female_lo = xr_open_lazy("%s/sims/ISMED-CLIM/V2511A/sims_model_V2511A_Portugal_female_lo.nc" %(DIR_DATA))
