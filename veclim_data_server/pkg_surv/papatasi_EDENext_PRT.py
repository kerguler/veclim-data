import json

from ..environ import DIR_DATA

print("Loading surv_papatasi_EDENext_PRT...",flush=True)

surv = json.load(open("%s/surveillance/EDENext/surv_papatasi_EDENext_PRT.json" %(DIR_DATA), "r"))