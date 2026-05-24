from ..environ import DIR_DATA
from ..functions import xr_open_lazy

print("Loading papatasi_V2511A_PRT...",flush=True)

import numpy
import xarray

from typing import Iterable, List, Sequence

import geopandas as gpd

from shapely.prepared import prep
from shapely.geometry import Point
from shapely import points

import datetime
import calendar

def mergeClose(up_times, down_times, sep):
    new_up_times = {}
    new_down_times = {}
    #
    for poly in up_times.keys():
        ups = list(up_times.get(poly, []))
        downs = list(down_times.get(poly, []))
        #
        if not ups:
            new_up_times[poly] = []
            new_down_times[poly] = []
            continue
        #
        merged_up = [ups[0]]
        merged_down = []
        cur_down = downs[0]
        #
        for u, d in zip(ups[1:], downs[1:]):
            if u - cur_down < sep:
                cur_down = max(cur_down, d)
            else:
                merged_down.append(cur_down)
                merged_up.append(u)
                cur_down = d
        #
        merged_down.append(cur_down)
        #
        # remove peaks with width < sep
        keep = [(u, d) for u, d in zip(merged_up, merged_down) if (d - u) >= sep]
        #
        new_up_times[poly] = [u for u, d in keep]
        new_down_times[poly] = [d for u, d in keep]
        #
    return new_up_times, new_down_times

def crossings_to_dict(cross_bool):
    # Stack poly/time -> a 1D MultiIndex and convert to a pandas Series (this computes)
    s = cross_bool.stack(evt=("poly", "time")).to_series()
    # Keep only True entries
    s = s[s]
    # Group by poly level and collect time values
    out = (
        s.index.to_frame(index=False)
         .groupby("poly", sort=False)["time"]
         .apply(list)
         .to_dict()
    )
    # Ensure all polys are present (even if no crossings)
    # (if you truly want every poly key)
    out = {str(poly): out.get(poly, []) for poly in cross_bool["poly"].values}
    #
    return out

def getCrossings(da,thresh=1.0,sep=0.0):
    above = da > thresh                       # (poly,time) bool
    prev  = above.shift(time=1)               # previous timestep (has NaN at first)
    #
    up_cross   = above & (prev == False)      # False -> True  (upward breach)
    down_cross = (~above) & (prev == True)    # True  -> False (downward breach)
    #
    up_times = crossings_to_dict(up_cross)
    down_times = crossings_to_dict(down_cross)
    #
    ids = sorted(set(up_times) | set(down_times))
    bad = [p for p in ids if len(up_times[p]) != len(down_times[p])]
    if bad:
        raise ValueError(f"Unpaired crossings for polygons: {bad[:10]}{'...' if len(bad)>10 else ''}")
    #
    if sep > 0.0:
        return mergeClose(up_times, down_times, sep)
    return up_times, down_times

def getPolyProp(polys,lon,lat):
    lonlat_crs = "EPSG:4326"
    pt = gpd.GeoSeries([Point(lon, lat)], crs=lonlat_crs).to_crs(polys.crs).iloc[0]
    #
    sidx = polys.sindex  # builds/uses R-tree (pygeos/shapely2 or rtree)
    cand_idx = list(sidx.intersection(pt.bounds))
    if not cand_idx:
        return None
    #
    cand = polys.iloc[cand_idx]
    #
    hits = cand[cand.geometry.contains(pt)]
    if hits.empty:
        return None
    if len(hits) > 1:
        raise ValueError(f"Multiple polygons matched ({len(hits)}).")
    #
    return hits.iloc[0]

def panelCube(base_nc, poly_ids):
    poly_ids = numpy.asarray(poly_ids)
    n_poly = poly_ids.size
    n_time = base_nc.sizes["time"]
    labels = list(base_nc.data_vars)
    #
    data_mat = numpy.full(
        (n_poly, n_time),
        numpy.nan,
        dtype=numpy.float32,
    )
    # Build output coords
    coords = {
        "poly": ("poly", poly_ids),
        "time": base_nc["time"].copy(),
    }
    # Copy scalar coords if present (same style as your example: 0-d coords)
    for c in ["step", "surface", "valid_time"]:
        if c in base_nc.coords:
            coords[c] = base_nc.coords[c].copy()
    #
    data_vars = {
        label: (("poly", "time"), data_mat)
        for label in labels
    }
    ds_out = xarray.Dataset(data_vars=data_vars, coords=coords)
    ds_out.attrs.update(base_nc.attrs)
    # Make conventions explicit if you want
    ds_out.attrs.setdefault("Conventions", base_nc.attrs.get("Conventions", "CF-1.7"))
    ds_out.attrs["zonal_aggregation"] = "polygon time-series means (poly,time) derived from original (y,x,time)"
    # Reasonable default compression/chunking for (poly,time)
    encoding = {
        label: {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
            # Chunk mostly along time for fast time-series access
            # Tune these depending on your n_poly / n_time
            "chunksizes": (min(n_poly, 512), min(n_time, 256)),
        }
        for label in labels
    }
    # ds_out.to_netcdf(out_path, encoding=encoding)
    return ds_out, encoding

def _in_interval(day: int, start: int, end: int, ndays: int = 365) -> bool:
    """
    Inclusive interval membership on a circular day-of-year axis.
    Assumes day/start/end are in [1, ndays].
    """
    if start <= end:
        return start <= day <= end
    # wrap-around case, e.g. 350..365 and 1..40
    return day >= start or day <= end

def _in_any_interval(day: int, ups: Sequence[int], downs: Sequence[int], ndays: int = 365) -> bool:
    if len(ups) != len(downs):
        raise ValueError("ups and downs must have the same length (paired intervals).")
    return any(_in_interval(day, u, d, ndays=ndays) for u, d in zip(ups, downs))

def classify_days(
    days: Iterable[int],
    up: Sequence[int],
    down: Sequence[int],
    peak_up: Sequence[int],
    peak_down: Sequence[int],
    ndays: int = 365,
    ) -> List[int]:
    """
    Return a list of {0,1,2} labels for each input day-of-year (1..365).
    Priority: peak (2) > active (1) > none (0).

    --- example usage ---
    days = [10, 40, 120, 200, 355]
    up = [30, 180]
    down = [150, 260]
    peak_up = [60, 210]
    peak_down = [90, 230]
    labels = classify_days(days, up, down, peak_up, peak_down)
    """
    out: List[int] = []
    for day in days:
        if not (1 <= day <= ndays):
            raise ValueError(f"Day {day} outside 1..{ndays}")
            #
        in_peak = _in_any_interval(day, peak_up, peak_down, ndays=ndays)
        if in_peak:
            out.append(2)
            continue
            #
        in_active = _in_any_interval(day, up, down, ndays=ndays)
        out.append(1 if in_active else 0)
        #
    return out

def times_to_dates(dict_times):
    year = datetime.date.today().year
    is_leap = calendar.isleap(year)
    out = {}
    for pid, times in dict_times.items():
        dates = []
        for doy in times:
            # shift by +1 day after Feb 28 in leap years
            if is_leap and doy >= 59:  # 59 = Feb 28 in 0-based index
                doy_shifted = doy + 1
            else:
                doy_shifted = doy
            date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=int(doy_shifted))
            dates.append(date)
        out[pid] = dates
    return out

def getWarning_v1(pid):
    ret = {
        'Season length (days)': (down_times[pid][-1]-up_times[pid][0]) if len(up_times[pid])>0 else numpy.nan,
        'Number of peaks': len(up_times[pid]),
        'Season start(s) (S0->S1)': "%s" %(", ".join(["%s" %(up_dates[pid][i].strftime("%d/%m/%Y")) for i in range(len(up_dates[pid]))])) if len(up_dates[pid])>0 else '',
        'Season end(s) (S1->S0)': "%s" %(", ".join(["%s" %(down_dates[pid][i].strftime("%d/%m/%Y")) for i in range(len(down_dates[pid]))])) if len(down_dates[pid])>0 else '',
        'Peak start(s) (S1->S2)': "%s" %(", ".join(["%s" %(peak_up_dates[pid][i].strftime("%d/%m/%Y")) for i in range(len(peak_up_dates[pid]))])) if len(peak_up_dates[pid])>0 else '',
        'Peak end(s) (S2->S1)': "%s" %(", ".join(["%s" %(peak_down_dates[pid][i].strftime("%d/%m/%Y")) for i in range(len(peak_down_dates[pid]))])) if len(peak_down_dates[pid])>0 else '',
        'Pre-season alert': "%s" %(", ".join(["%s" %((up_dates[pid][i] + datetime.timedelta(days=int(-14))).strftime("%d/%m/%Y")) for i in [0]])) if len(up_dates[pid])>0 else '',
        'Start-of-season alert': "%s" %(", ".join(["%s" %((peak_up_dates[pid][i] + datetime.timedelta(days=int(-14))).strftime("%d/%m/%Y")) for i in [0]])) if len(peak_up_dates[pid])>0 else '',
        'Low activity alert': "%s" %(", ".join(["%s" %((peak_down_dates[pid][i] + datetime.timedelta(days=int(0))).strftime("%d/%m/%Y")) for i in [-1]])) if len(peak_down_dates[pid])>0 else '',
        'End-of-season alert': "%s" %(", ".join(["%s" %((down_dates[pid][i] + datetime.timedelta(days=int(14))).strftime("%d/%m/%Y")) for i in [-1]])) if len(down_dates[pid])>0 else '',
        'Legend': """
<div class="legend">
  <h3>Legend</h3>

  <p>
    This output summarises the predicted seasonal dynamics of vector activity at a given location
    (longitude, latitude, area code, and area name). The season length (in days) denotes the total
    duration of the active period, while the number of peaks indicates how many high-activity
    periods occur within the season.
  </p>

  <p><strong>Seasonal phases are defined using discrete risk levels:</strong></p>
  <ul>
    <li><strong>S0:</strong> low risk</li>
    <li><strong>S1:</strong> moderate risk</li>
    <li><strong>S2:</strong> high risk</li>
  </ul>

  <p><strong>Transitions between these states define key dates:</strong></p>
  <ul>
    <li><strong>Season start (S0 &rarr; S1):</strong> onset of the active season</li>
    <li><strong>Season end (S1 &rarr; S0):</strong> termination of the active season</li>
    <li><strong>Peak start (S1 &rarr; S2):</strong> beginning of high-risk period</li>
    <li><strong>Peak end (S2 &rarr; S1):</strong> return from high to moderate risk</li>
  </ul>

  <p>
    Each transition is accompanied by notification dates designed to provide staged early warnings
    of seasonal changes: a pre-season alert is issued two weeks before the first season start
    (S0 &rarr; S1), followed by a start-of-season alert two weeks before the peak season onset
    (if applicable). A low-activity alert is issued at the end of the peak period
    (S2 &rarr; S1), if applicable, and an end-of-season alert is given two weeks after the final
    transition back to no-risk conditions (S1 &rarr; S0).
  </p>
</div>
"""
    }
    return ret

def get_warning(pid):
    DATE_START             = up_dates[pid][0] if len(up_dates[pid])>0 else None
    DATE_END               = down_dates[pid][-1] if len(down_dates[pid])>0 else None
    DATE_PEAK_START        = peak_up_dates[pid][0] if len(peak_up_dates[pid])>0 else None
    DATE_PEAK_END          = peak_down_dates[pid][-1] if len(peak_down_dates[pid])>0 else None
    DATE_LOW_START         = DATE_PEAK_END
    DATE_LOW_END           = DATE_END
    #
    PRE_SEASON_ALERT       = DATE_START + datetime.timedelta(days=int(-14)) if DATE_START is not None else None
    START_OF_SEASON_ALERT  = DATE_START
    PEAK_ACTIVITY_REMINDER = DATE_PEAK_START
    LOW_ACTIVITY_REMINDER  = DATE_LOW_START
    END_OF_SEASON_ALERT    = DATE_END
    #
    ret = {
        'Area code': pid,
        'Dates': {
            'DATE_START': DATE_START.strftime("%d/%m/%Y") if DATE_START is not None else '',
            'DATE_END': DATE_END.strftime("%d/%m/%Y") if DATE_END is not None else '',
            'DATE_PEAK_START': DATE_PEAK_START.strftime("%d/%m/%Y") if DATE_PEAK_START is not None else '',
            'DATE_PEAK_END': DATE_PEAK_END.strftime("%d/%m/%Y") if DATE_PEAK_END is not None else '',
            'DATE_LOW_START': DATE_LOW_START.strftime("%d/%m/%Y") if DATE_LOW_START is not None else '',
            'DATE_LOW_END': DATE_LOW_END.strftime("%d/%m/%Y") if DATE_LOW_END is not None else '',
        },
        'Notifications': {}
    }
    if (PRE_SEASON_ALERT is not None) and (DATE_START is not None) and (DATE_END is not None):
        ret['Notifications']['PRE_SEASON_ALERT'] = {
            'Send': PRE_SEASON_ALERT.strftime("%d/%m/%Y"),
            'Subject': "Sand fly season starts soon. Protect your dog!",
            'Message': "Our Early Warning and Response System predicts that sand flies will become active in your area on %s." %(DATE_START.strftime("%d/%m/%Y"))
        }
    if (START_OF_SEASON_ALERT is not None) and (DATE_START is not None) and (DATE_END is not None):
        ret['Notifications']['START_OF_SEASON_ALERT'] = {
            'Send': START_OF_SEASON_ALERT.strftime("%d/%m/%Y"),
            'Subject': "Sand fly season is now active. Act today!",
            'Message': "Sand flies are now active in your area (%s, %s)." %(DATE_START.strftime("%d/%m/%Y"), DATE_END.strftime("%d/%m/%Y"))
        }
    if (PEAK_ACTIVITY_REMINDER is not None) and (DATE_PEAK_START is not None) and (DATE_PEAK_END is not None):
        ret['Notifications']['PEAK_ACTIVITY_REMINDER'] = {
            'Send': PEAK_ACTIVITY_REMINDER.strftime("%d/%m/%Y"),
            'Subject': "Peak sand fly activity expected between %s and %s." %(DATE_PEAK_START.strftime("%d/%m/%Y"), DATE_PEAK_END.strftime("%d/%m/%Y")),
            'Message': "This is the highest-risk period for transmission. Please strictly follow all preventive measures."
        }
    if (LOW_ACTIVITY_REMINDER is not None) and (DATE_LOW_START is not None) and (DATE_LOW_END is not None):
        ret['Notifications']['LOW_ACTIVITY_REMINDER'] = {
            'Send': LOW_ACTIVITY_REMINDER.strftime("%d/%m/%Y"),
            'Subject': "Reduced sand fly activity expected between %s and %s." %(DATE_LOW_START.strftime("%d/%m/%Y"), DATE_LOW_END.strftime("%d/%m/%Y")),
            'Message': "Sand fly activity is expected to be low during this period. However, transmission is still possible. Please continue using repellents and keeping your dog indoors at night."
        }
    if (END_OF_SEASON_ALERT is not None) and (DATE_END is not None):
        ret['Notifications']['END_OF_SEASON_ALERT'] = {
            'Send': END_OF_SEASON_ALERT.strftime("%d/%m/%Y"),
            'Subject': "Sand fly season is over. Time for post-season testing!",
            'Message': "Sand fly infection risk is expected to decrease by %s. Please book a post-season Leishmaniasis test if your veterinarian recommends it. Thank you for protecting animal and public health." %(DATE_END.strftime("%d/%m/%Y"))
        }
    #
    return ret

class dbPortugal:
    def __init__(self,var,filename="",nc=None,shapefile="",verbose=False):
        self.var = var
        self.filename = filename        
        self.nc = nc
        self.shapefile = shapefile
        if (self.filename == "") and (self.nc is None):
                raise ValueError(f"Please supply either the simulation matrix or the panelCube file.")
        #
        self.polys = gpd.read_file(self.shapefile)
        self.polys = self.polys.to_crs("EPSG:4326")
        self.poly_ids = self.polys['Official_Co']
        #
        if self.filename == "":
            self.mat, self.encoding = panelCube(self.nc,self.poly_ids)
            means, grids = self.setmat(verbose=verbose)
        else:
            self.mat = xr_open_lazy(self.filename)
            self.encoding = self.mat[var].encoding
            means = self.mat[var].mean("time", skipna=True)
        #
        self.polys['Means'] = means
        #
    def setmat(self,verbose=False):
        lon2d = (((self.nc["longitude"].values + 180) % 360) - 180)
        lat2d = self.nc["latitude"].values
        #
        grids = {}
        means = []
        for i in range(len(self.polys)):
            if verbose:
                print("Processing %d of %d..." %(i+1,len(self.polys)))
            #
            poly = self.polys.geometry.iloc[i]
            label = self.polys['Official_Co'].iloc[i,0]
            P = prep(poly)
            #
            # fast bbox filter first
            minx, miny, maxx, maxy = poly.bounds
            cand = (lon2d >= minx) & (lon2d <= maxx) & (lat2d >= miny) & (lat2d <= maxy)
            #
            if not numpy.any(cand):
                means.append(numpy.nan)
                continue
            #
            yy, xx = numpy.where(cand)
            pts = points(lon2d[yy, xx], lat2d[yy, xx])      # vectorized
            #
            inside_cand = numpy.array([P.contains(p) for p in pts])
            #
            yx = numpy.column_stack([yy[inside_cand], xx[inside_cand]])
            grids[label] = yx
            #
            tmp = [
                self.nc[self.var].isel(y=y, x=x).load().values 
                for y,x in yx
                ]
            if len(tmp) == 0 or numpy.all(numpy.isnan(tmp)):
                means.append(numpy.nan)
            else:
                means.append(numpy.nanmean(tmp))
                self.mat[self.var].loc[dict(poly=label)] = numpy.nanmean(tmp,axis=0)
            #
        return means, grids

db = dbPortugal("newegg",
                filename="%s/sims/ISMED-CLIM/V2511A_PRT/sims_model_V2511A_Portugal_newegg_mean_poly.nc" %DIR_DATA,
                shapefile="%s/sims/ISMED-CLIM/V2511A_PRT/georef-portugal-concelho-millesime.shp" %DIR_DATA,
                verbose=False)

up_times, down_times = getCrossings(db.mat["newegg"], 
                                    thresh=1.0, 
                                    sep=14.0)
up_dates = times_to_dates(up_times)
down_dates = times_to_dates(down_times)

peak_up_times, peak_down_times = getCrossings(db.mat["newegg"], 
                                              thresh=10000.0, 
                                              sep=14.0)
peak_up_dates = times_to_dates(peak_up_times)
peak_down_dates = times_to_dates(peak_down_times)