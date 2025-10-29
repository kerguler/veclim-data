import numpy

from ..environ import DIR_DATA
from ..functions import xr_open_lazy, cache_npy

print("Loading NASA...",flush=True)

models = ['ACCESS-CM2','ACCESS-ESM1-5','EC-Earth3']
ssps = ['ssp245','ssp585']
decade = "2090-2099"
dlabel = "2090_to_2099"

def calc_colegg(ssp):
    return numpy.nanmean([
        xr_open_lazy("%s/NASA_decadal_%s_%s/%s/sims_vector08c_Q4.a100+1_NASA_decadal_%s_%s_%s_colegg.nc" %(dr,model,ssp,dlabel,model,ssp,dlabel))['colegg']
        for model in models
    ], axis=0)

def calc_iouts(ssp):
    return numpy.nanmean([
        xr_open_lazy("%s/NASA_decadal_%s_%s/p4000r100w60/%s/sims_vector08c_Q4.a100+1_chikv_QI_NASA_decadal_%s_%s_p4000r100w60_%s_iouts.nc" %(dr,model,ssp,dlabel,model,ssp,dlabel))['iouts']
        for model in models
    ], axis=0)

def calc_pouts(ssp):
    return numpy.nanmean([
        xr_open_lazy("%s/NASA_decadal_%s_%s/p4000r100w60/%s/sims_vector08c_Q4.a100+1_chikv_QI_NASA_decadal_%s_%s_p4000r100w60_%s_pouts.nc" %(dr,model,ssp,dlabel,model,ssp,dlabel))['pouts']
        for model in models
    ], axis=0)

dr = "%s/sims/vector08c_Q4.a100+1" %(DIR_DATA)
colegg = {
    ssp: cache_npy("annualNASA_colegg_%s.npy" %ssp,calc_colegg,ssp)
    for ssp in ssps
}

dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI" %(DIR_DATA)
iouts = {
    ssp: cache_npy("annualNASA_iouts_%s.npy" %ssp,calc_iouts,ssp)
    for ssp in ssps
}
pouts = {
    ssp: cache_npy("annualNASA_pouts_%s.npy" %ssp,calc_pouts,ssp)
    for ssp in ssps
}

popsize, reps, win = [4000, 100, 60]