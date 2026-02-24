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

def mergeClose(up_times, down_times, sep):
    cmp = lambda gap: gap < sep
    #
    new_up_times = {}
    new_down_times = {}
    #
    ids = list(up_times.keys())
    for poly in ids:
        ups = list(up_times.get(poly, []))
        downs = list(down_times.get(poly, []))
        #
        if not ups:
            new_up_times[poly] = []
            new_down_times[poly] = []
            continue
        #
        keep_up = [ups[0]]
        keep_down = []
        #
        cur_down = downs[0]
        #
        for u, d in zip(ups[1:], downs[1:]):
            gap = u - cur_down
            #
            if cmp(gap):
                # merge: drop this up, extend current down
                cur_down = max(cur_down, d)
            else:
                # close current peak, start new one
                keep_down.append(cur_down)
                keep_up.append(u)
                cur_down = d
        #
        keep_down.append(cur_down)
        #
        new_up_times[poly] = keep_up
        new_down_times[poly] = keep_down
        #
    return new_up_times, new_down_times

def crossings_to_dict(cross_bool):
    stacked = cross_bool.stack(evt=("poly", "time"))
    stacked = stacked.where(stacked, drop=True)
    mi = stacked.indexes["evt"]  # MultiIndex (poly, time)
    out = {}
    for poly in mi.levels[0]:
        # select times for this poly (might be empty)
        times = mi[mi.get_level_values("poly") == poly].get_level_values("time").to_list()
        out[poly] = times
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

        in_peak = _in_any_interval(day, peak_up, peak_down, ndays=ndays)
        if in_peak:
            out.append(2)
            continue

        in_active = _in_any_interval(day, up, down, ndays=ndays)
        out.append(1 if in_active else 0)

    return out

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
        self.poly_ids = self.polys['Official_Co'].iloc[:,0]
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
                filename="%s/sims/ISMED-CLIM/V2511A_PRT/sims_model_V2511A_Portugal_newegg_poly.nc" %DIR_DATA,
                shapefile="%s/sims/ISMED-CLIM/V2511A_PRT/georef-portugal-concelho-millesime.shp" %DIR_DATA,
                verbose=False)

up_times, down_times = getCrossings(db.mat["newegg"], 
                                    thresh=1.0, 
                                    sep=14.0)

peak_up_times, peak_down_times = getCrossings(db.mat["newegg"], 
                                              thresh=1000.0, 
                                              sep=14.0)