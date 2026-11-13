import json

from ..environ import DIR_DATA

print("Loading surv_sandfly_EDENext...",flush=True)

survPRT = json.load(open("%s/surveillance/EDENext/surv_papatasi_EDENext_PRT.json" %(DIR_DATA), "r"))

surv = json.load(open("%s/surveillance/EDENext/surv_EDENext_CERRA.json" %(DIR_DATA), "r"))
