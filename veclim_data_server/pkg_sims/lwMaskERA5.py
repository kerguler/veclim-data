import numpy
import json

from ..environ import DIR_DATA

lwmask = numpy.load("%s/clim/lwmask_0.1lw.npy" %(DIR_DATA),mmap_mode='r')
lonlat = json.load(open("%s/clim/lonlat.json" %(DIR_DATA),"r"))
latitude = numpy.array(lonlat['lat'])
longitude = numpy.array(lonlat['lon'])

def island(loni,lati):
    return ~lwmask[lati,loni]
