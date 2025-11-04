import os
import numpy
import pandas
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

def isel(nc,lon=None,lat=None):
    kw = {}
    if type(lon) != type(None):
        kw['lon'] = lon
    if type(lat) != type(None):
        kw['lat'] = lat
    return nc.isel(**kw).load().values

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

def calc_cut(vec,lim,lab):
    if hasattr(vec, '__iter__'):
        return pandas.cut(vec, bins=lim, include_lowest=False, right=False, labels=lab).tolist()
    return lab[numpy.where(numpy.array(lim) > vec)[0][0] - 1]

def remove_feb29(vec,days,isFeb29,fill=0.0):
    tmp = vec[days-1]
    if len(isFeb29) == 0:
        if fill == None:
            return [None if numpy.isnan(a) else a for a in tmp]
        else:
            return numpy.nan_to_num(tmp,nan=fill).tolist()
    tmp[isFeb29+1] = 0.5*(tmp[isFeb29]+tmp[isFeb29+1])
    if fill == None:
        return [None if numpy.isnan(a) else a for a in numpy.delete(tmp,isFeb29)]
    else:
        return numpy.nan_to_num(numpy.delete(tmp,isFeb29),nan=fill).tolist()

from datetime import datetime, timedelta
Feb29 = datetime(2020,2,29).timetuple().tm_yday

not_clim_keys = [
    "inv",
    "lon",
    "lat",
    "day0",
    "day1",
]

def calc_list(clms):
    return {
        key: numpy.hstack([clm[key] for clm in clms]).tolist()
        for key in clms[0] if not key in not_clim_keys
    }

def calc_mean(clms):
    return {
        key: float(numpy.nanmean(numpy.hstack([clm[key] for clm in clms])))
        for key in clms[0] if not key in not_clim_keys
    }

def get_clim(get_days, loni, lati, pr0, pr1, ts=False):
    if ts:
        calc_fun = calc_list
    else:
        calc_fun = calc_mean
    #
    ret = calc_fun([get_days(loni,lati,pr0,pr1)])
    #
    return ret

def get_dates(date0, date1=False, ts=False):
    valid = True
    if (not date1) or (date1 < date0):
        date1 = date0
        valid = False
    #
    isFeb29 = []
    idates = []
    ddates = []
    while date0 <= date1:
        ddates.append(date0)
        idates.append(date0.timetuple().tm_yday-1 if is_leap_year(date0.year) and date0.timetuple().tm_yday>Feb29 else date0.timetuple().tm_yday)
        if date0.month==2 and date0.day==29:
            isFeb29.append(len(ddates)-1)
            ddates.append(date0)
            idates.append(Feb29-1)
        date0 += timedelta(days=1)
    #
    ddates = numpy.array(ddates)
    idates = numpy.array(idates)
    isFeb29 = numpy.array(isFeb29)
    #
    return {
        "dates": ddates,
        "days": idates,
        "isFeb29": isFeb29,
        "date0": ddates[0].strftime("%Y-%m-%d"),
        "date1": ddates[-1].strftime("%Y-%m-%d"),
        "valid": int(valid)
    }