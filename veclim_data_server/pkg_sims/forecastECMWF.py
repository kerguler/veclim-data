import os
import sys
import re

from ..environ import DIR_DATA
from ..functions import xr_open_lazy, stdwarn

dr = "%s/sims/vector08c_Q4.a100+1/ECMWF" %(DIR_DATA)

fld = sorted(os.listdir(dr))[-1]
overlap = [re.sub(r"\_","-",s) for s in re.findall(r"\d\d\d\d\_\d+\_\d+",fld)]

if (len(overlap)!=2):
    stdwarn("Vector forecast directory is misworded! %s" %fld)
    sys.exit(1)

colegg = xr_open_lazy("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_colegg.nc" %(dr,fld,fld))['colegg']
colK = xr_open_lazy("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_colK.nc" %(dr,fld,fld))['colK']
coln2 = xr_open_lazy("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_coln2.nc" %(dr,fld,fld))['coln2']
coln4f = xr_open_lazy("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_coln4f.nc" %(dr,fld,fld))
coln4f = coln4f['coln4f']
dates = coln4f['time']

dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI/ECMWF" %(DIR_DATA)

fld = sorted(os.listdir(dr))[-1]
overlap = [re.sub(r"\_","-",s) for s in re.findall(r"\d\d\d\d\_\d+\_\d+",fld)]

print("Reading forecast...", overlap)

if (len(overlap)!=2):
    stdwarn("Disease forecast directory is misworded! %s" %fld)
    sys.exit(1)

if (overlap[0] != overlap[0]) or (overlap[1] != overlap[1]):
    stdwarn("Dates of vector and disease forcasts mismatch! %s != %s or %s != %s" %(overlap[0],overlap[0],overlap[1],overlap[1]))
    sys.exit(1)

tmp = re.findall(r"p(\d+?)r(\d+?)w(\d+?)",fld)
if len(tmp)!=1:
    stdwarn("Disease forecast configuration is missing or misworded! %s" %fld)
    sys.exit(1)

popsize, reps, win = tmp[0]

pouts = xr_open_lazy("%s/%s/sims_vector08c_Q4.a100+1_chikv_QI_ECMWF_%s_pouts.nc" %(dr,fld,fld))['pouts']
iouts = xr_open_lazy("%s/%s/sims_vector08c_Q4.a100+1_chikv_QI_ECMWF_%s_iouts.nc" %(dr,fld,fld))['iouts']
idates = iouts['time']
iouts = iouts

def getOverlap():
    return overlap

def getTransConf():
    return {
        'popsize': popsize, 
        'reps': reps, 
        'win': win
    }
