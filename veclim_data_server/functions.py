import os
import numpy
import pandas
import xarray

from veclim_data_server.environ import DIR_CACHE

from datetime import datetime, timedelta
Feb29 = datetime(2020,2,29).timetuple().tm_yday

not_clim_keys = [
    "inv",
    "lon",
    "lat",
    "day0",
    "day1",
]

def cache_npy(filename, func, *args, **kwargs):
    fildir = "%s/%s" %(DIR_CACHE,filename)
    if os.path.exists(fildir):
        try:
            mat = numpy.load(fildir,mmap_mode='r')
            return mat
        except:
            pass
    try:
        mat = func(*args, **kwargs)
        numpy.save(fildir, mat)
    except Exception as e:
        print(e, flush=True)
        print("ERROR: Failed to create %s" %filename, flush=True)
        return []
    return mat

def cache_ncdf(filename, func, *args, **kwargs):
    fildir = "%s/%s" %(DIR_CACHE,filename)
    if os.path.exists(fildir):
        try:
            mat = xr_open_lazy(fildir)
            return mat
        except:
            pass
    try:
        mat = func(*args, **kwargs)
        mat.to_netcdf(fildir,
                engine="netcdf4", 
                compute=True,
                encoding={
                    key: {
                        "zlib": True, 
                        "complevel": 9
                    } for key in mat.data_vars
                })
    except Exception as e:
        print(e, flush=True)
        print("ERROR: Failed to create %s" %filename, flush=True)
        return []
    return mat

def xr_open_lazy(path, engine=None, chunks="auto"):
    kw = {"chunks": chunks}
    if engine:
        kw["engine"] = engine
    return xarray.open_dataset(path, **kw)

def isel(nc, lon=None, lat=None, x=None, y=None, time=None,
         ret_df=False, loni=None, lati=None):

    kw_isel = {}
    kw_sel = {}

    # Detect coordinate/dimension names
    lon_name = next(
        (k for k in ('lon', 'longitude') if k in nc.dims), None
    )
    lat_name = next(
        (k for k in ('lat', 'latitude') if k in nc.dims), None
    )

    if loni is not None:
        kw_isel[lon_name] = loni
    elif lon is not None:
        kw_sel[lon_name] = lon

    if lati is not None:
        kw_isel[lat_name] = lati
    elif lat is not None:
        kw_sel[lat_name] = lat

    for key, val in [('x', x), ('y', y), ('time', time)]:
        if val is not None:
            kw_isel[key] = val

    if kw_sel:
        nc = nc.sel(**kw_sel, method='nearest')

    tmp = nc.isel(**kw_isel).load()

    return tmp if ret_df else tmp.values

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

def remove_feb29(vec, days, isFeb29, fill=0.0, mean=True, tolist=True):
    vec = numpy.asarray(vec, dtype=float)
    days = numpy.asarray(days, dtype=int)
    isFeb29 = numpy.asarray(isFeb29, dtype=int)
    # Build output by indexing vec with day indices (days is 0-based)
    tmp = numpy.take(vec, indices=days, axis=-1)
    #
    if isFeb29.size:
        if mean:
            # For Feb29 positions, days[pos] is the collapsed Feb28 index
            feb28_idx = days[isFeb29]
            mar01_idx = numpy.clip(feb28_idx + 1, 0, vec.shape[-1] - 1)
            #
            feb28_vals = numpy.take(vec, indices=feb28_idx, axis=-1)
            mar01_vals = numpy.take(vec, indices=mar01_idx, axis=-1)
            #
            tmp[..., isFeb29] = 0.5 * (feb28_vals + mar01_vals)
            #
    if fill is not None:
        rv = numpy.nan_to_num(tmp, nan=fill)
        return rv.tolist() if tolist else rv
        #
    flat = tmp.reshape(-1)
    flat = pandas.Series(flat).where(~numpy.isnan(flat), None).to_numpy()
    rv = flat.reshape(tmp.shape)
    return rv.tolist() if tolist else rv


def isel_feb29(nc, idates=[], isFeb29=[], lon=None, lat=None, x=None, y=None,
               fill=0.0, mean=True, tolist=True, loni=None, lati=None):
    #
    days = numpy.asarray(idates, dtype=int)
    isFeb29 = numpy.asarray(isFeb29, dtype=int)
    #
    # Preserve original order, removing duplicates
    indices = list(dict.fromkeys(days))
    #
    # Include additional days needed for interpolation
    if mean and isFeb29.size:
        indices = list(dict.fromkeys(
            indices + list(days[isFeb29] + 1)
        ))
    #
    # Load selected time steps
    vec = isel(
        nc, lon=lon, lat=lat, x=x, y=y,
        loni=loni, lati=lati, time=indices,
        ret_df=True
    )
    #
    # Convert Dataset to DataArray if necessary
    if isinstance(vec, xarray.Dataset):
        if len(vec.data_vars) != 1:
            raise ValueError(
                "Expected a Dataset containing one data variable."
            )
        vec = next(iter(vec.data_vars.values()))
    #
    # Ensure time is the last dimension
    vec = vec.transpose(
        *[d for d in vec.dims if d != 'time'],
        'time'
    ).values
    #
    # Map original indices to selected positions
    lookup = {idx: pos for pos, idx in enumerate(indices)}
    positions = [lookup[d] for d in days]
    #
    tmp = numpy.take(vec, positions, axis=-1)
    #
    # Interpolate February 29
    if mean and isFeb29.size:
        next_positions = [lookup[d + 1] for d in days[isFeb29]]
        #
        tmp[..., isFeb29] = 0.5 * (
            tmp[..., isFeb29] +
            numpy.take(vec, next_positions, axis=-1)
        )
    #
    if fill is not None:
        tmp = numpy.nan_to_num(tmp, nan=fill)
        return tmp.tolist() if tolist else tmp
    #
    flat = tmp.reshape(-1)
    flat = pandas.Series(flat).where(~numpy.isnan(flat), None).to_numpy()
    rv = flat.reshape(tmp.shape)
    #
    return rv.tolist() if tolist else rv

def calc_list(clms):
    return {
        key: numpy.hstack([clm[key] for clm in clms]).tolist()
        for key in clms[0] if not key in not_clim_keys
    }

def calc_mean(clms):
    return {
        key: float(numpy.nanmean(
            numpy.hstack([clm[key] for clm in clms])
            ))
        for key in clms[0] if not key in not_clim_keys
    }

def get_clim(get_days, loni, lati, pr0, pr1, ts=False):
    if ts:
        calc_fun = calc_list
    else:
        calc_fun = calc_mean
    #
    ret = calc_fun([get_days(loni, lati, pr0, pr1)])
    #
    return ret

def get_dates(date0, date1=None, ts=True):
    """
    Datetime-based date generator.

    Returns
    -------
    dict with:
        dates : numpy array of datetime.datetime
        day_index_365 : numpy int array (0..364)
        isFeb29 : numpy int array (positions in date list)
        date0 : ISO string
        date1 : ISO string
        valid : int
    """
    date0 = pandas.Timestamp(date0)
    date1 = pandas.Timestamp(date1) if date1 is not None else None
    #
    valid = 1
    if date1 is None or date1 < date0:
        date1 = date0
        valid = 0
    # Internal pandas index
    dates_index = pandas.date_range(start=date0, end=date1, freq="D")
    # 1..366
    dayofyear = dates_index.dayofyear.to_numpy()
    # Leap year mask (already numpy array)
    is_leap = dates_index.is_leap_year
    # Collapse leap years to 365-day calendar
    day_365 = dayofyear.copy()
    day_365[is_leap & (dayofyear >= 60)] -= 1
    days = day_365 - 1  # convert to 0-based index
    #
    isFeb29 = numpy.flatnonzero((dates_index.month == 2) & (dates_index.day == 29))
    # Convert to datetime.datetime objects (NOT numpy datetime64)
    dates = dates_index.to_pydatetime()
    #
    return {
        "dates": numpy.array(dates, dtype=object),
        "days": days.astype(numpy.int32),
        "isFeb29": isFeb29.astype(numpy.int32),
        "date0": dates[0].strftime("%Y-%m-%d"),
        "date1": dates[-1].strftime("%Y-%m-%d"),
        "valid": int(valid),
    }
