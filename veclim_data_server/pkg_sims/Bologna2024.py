from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading Bologna2024...",flush=True)

models = [
    "albopictus (sPop)",
    "ArboCartoR",
    "dynamAedes",
    "Metelmann et al. (2019)",
    "VECTRI",
    "AedesDDE",
    "Stacked Machine Learning"
]

ncs = [
    xr_open_lazy("%s/sims/Bologna2024/sims_Bologna2024_%d_%s.nc" %(DIR_DATA,id,model))
    for id, model in enumerate(models)
]

lons = ncs[0].lon.load().values
lats = ncs[0].lat.load().values
dates = ncs[0].yrwk.load().values