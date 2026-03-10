import json

from ..environ import DIR_DATA

print("Loading surv_papatasi_V2511A_PRT...",flush=True)

surv = json.load(open("%s/surveillance/ISMED-CLIM/surv_papatasi_V2511A_PRT.json" %(DIR_DATA), "r"))