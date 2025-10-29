import numpy

from ..environ import DIR_DATA
from ..functions import xr_open_lazy

decade = "2010-2019"
dlabel = "2010_to_2019"
dr = "%s/sims/vector08c_Q4.a100+1/ERA5_single_levels_decadal/%s" %(DIR_DATA,dlabel)
colegg = xr_open_lazy("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_colegg.nc" %(dr,dlabel))['colegg']
colK = xr_open_lazy("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_colK.nc" %(dr,dlabel))['colK']
coln2 = xr_open_lazy("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_coln2.nc" %(dr,dlabel))['coln2']
coln4f = xr_open_lazy("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_coln4f.nc" %(dr,dlabel))['coln4f']

dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI/ERA5_single_levels_decadal/p4000r100w60/%s" %(DIR_DATA,dlabel)
iouts = xr_open_lazy("%s/sims_vector08c_Q4.a100+1_chikv_QI_ERA5_single_levels_decadal_p4000r100w60_%s_iouts.nc" %(dr,dlabel))['iouts']
pouts = xr_open_lazy("%s/sims_vector08c_Q4.a100+1_chikv_QI_ERA5_single_levels_decadal_p4000r100w60_%s_pouts.nc" %(dr,dlabel))['pouts']

tmp = numpy.genfromtxt("%s/surveillance/Italy2008/coord_albopictus_Italy2008.csv" %(DIR_DATA), delimiter=',',names=True)
acc = {}
for row in tmp:
    loni = int(row['loni'])
    lati = int(row['lati'])
    if not lati in acc:
        acc[lati] = {}
    acc[lati][loni] = 1

def accuracy(lati,loni):
    if not lati in acc or not loni in acc[lati]:
        return 0
    return acc[lati][loni]