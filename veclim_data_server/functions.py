import os
import numpy
import xarray

def cache_npy(filename, func, *args, **kwargs):
    if os.path.exists(filename):
        try:
            mat = numpy.load(filename,mmap_mode='r')
            return mat
        except:
            pass
    print("Caching %s" %filename, flush=True)
    try:
        mat = func(*args, **kwargs)
        numpy.save(filename, mat)
    except:
        print("Failed to create %s" %filename, flush=True)
        return []
    return mat

def xr_open_lazy(path, engine=None, chunks="auto"):
    kw = {"chunks": chunks}
    if engine:
        kw["engine"] = engine
    return xarray.open_dataset(path, **kw)

def stdwarn(msg):
    print("ERROR",msg,flush=True)

def is_leap_year(y):
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False

def getIndex(vec,v):
    return numpy.argmin(numpy.abs(vec-v))

def daylength(lat,day):
    """
    Translated from the daylength function of the geosphere package of R
    lat: latitude in degree decimal (float)
    day: datetime.date or day of the year (integer)
    """
    from datetime import date
    if isinstance(day,date):
        day = day.timetuple().tm_yday
    pi180 = numpy.pi / 180.0
    P = numpy.arcsin(0.39795 * numpy.cos(0.2163108 + 2 * numpy.arctan(0.9671396 * numpy.tan(0.0086 * (day - 186)))))
    a = (numpy.sin(0.8333 * pi180) + numpy.sin(lat * pi180) * numpy.sin(P))/(numpy.cos(lat * pi180) * numpy.cos(P))
    a = numpy.min([numpy.max([a, -1]), 1])
    DL = 24 - (24.0/numpy.pi) * numpy.arccos(a)
    return(DL)