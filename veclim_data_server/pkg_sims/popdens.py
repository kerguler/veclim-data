import numpy

from ..environ import DIR_DATA

print("Loading popdens...",flush=True)

pop = numpy.load("%s/clim/SEDAC/gpw_v4_population_density_adjusted_fromHiRes_0.1lwmask_2010_to_2020.npy" %(DIR_DATA),mmap_mode='r')