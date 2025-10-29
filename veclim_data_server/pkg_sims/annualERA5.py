import numpy

from ..environ import DIR_DATA
from ..functions import xr_open_lazy, daylength

print("Loading annual ERA5...",flush=True)

decade = "2010-2019"
dr = "%s/clim/ERA5/ERA5_single_levels_decadal/2010_to_2019" %(DIR_DATA)
atemp = xr_open_lazy("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_temperature.nc" %(dr))
tf = 365
latitude = atemp['lat']
photo = numpy.array([[daylength(lat, d) for d in numpy.arange(tf)] for lat in latitude])
atemp = atemp['2m_temperature']
atmin = xr_open_lazy("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_temperature_min.nc" %(dr))['2m_temperature_min']
atmax = xr_open_lazy("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_temperature_max.nc" %(dr))['2m_temperature_max']
rehum = xr_open_lazy("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_relative_humidity.nc" %(dr))['2m_relative_humidity']
precp = xr_open_lazy("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_total_precipitation.nc" %(dr))['total_precipitation']
soilw = xr_open_lazy("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_volumetric_soil_water_layer_1.nc" %(dr))['volumetric_soil_water_layer_1']