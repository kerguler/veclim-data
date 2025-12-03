from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading surv_Bologna2024...",flush=True)

surv = xr_open_lazy("%s/Bologna2024/surv_Bologna2024.nc" %(DIR_DATA))