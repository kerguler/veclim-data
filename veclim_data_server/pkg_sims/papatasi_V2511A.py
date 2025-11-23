import numpy
from scipy.spatial import cKDTree

from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading papatasi_V2511A...",flush=True)

class gridCERRA:
    def __init__(self, filename=""):
        filename = filename if filename else "%s/clim/CERRA/clim_CERRA_lsmask.nc" %(DIR_DATA)
        self.ds = xr_open_lazy(filename)
        #
        self.lon = self.ds['longitude'].values
        # Convert from (0,360) to (-180,180)
        self.lon = ((self.lon + 180) % 360) - 180
        #
        self.lat = self.ds['latitude'].values
        self.points = numpy.column_stack((self.lon.ravel(), self.lat.ravel()))
        self.tree = cKDTree(self.points)
        #
    def geti(self,lon,lat):
        _, idx = self.tree.query([lon, lat])
        y, x = numpy.unravel_index(idx, self.lon.shape)
        return {
            'x': int(x),
            'y': int(y),
            'lon': float(self.lon[y,x]),
            'lat': float(self.lat[y,x])
        }
        #
    def get(self,lon,lat):
        ii = self.geti(lon,lat)
        ret = {
            'lon': numpy.float64(self.lon[ii['y'],ii['x']]),
            'lat': numpy.float64(self.lat[ii['y'],ii['x']])
        }
        return ret

female_md = xr_open_lazy("%s/sims/ISMED-CLIM/V2511A/sims_model_V2511A_Portugal_female_md.nc" %(DIR_DATA))
female_hi = xr_open_lazy("%s/sims/ISMED-CLIM/V2511A/sims_model_V2511A_Portugal_female_hi.nc" %(DIR_DATA))
female_lo = xr_open_lazy("%s/sims/ISMED-CLIM/V2511A/sims_model_V2511A_Portugal_female_lo.nc" %(DIR_DATA))

grid = gridCERRA()