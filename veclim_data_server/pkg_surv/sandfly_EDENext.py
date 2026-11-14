import json
import numpy

import geopandas as gpd
from shapely.geometry import box

from ..functions import getIndex
from veclim_data_server import pkg_sims

from ..environ import DIR_DATA

print("Loading surv_sandfly_EDENext...",flush=True)

survPRT = json.load(open("%s/surveillance/EDENext/surv_papatasi_EDENext_PRT.json" %(DIR_DATA), "r"))

surv = json.load(open("%s/surveillance/EDENext/surv_EDENext_CERRA.json" %(DIR_DATA), "r"))

lonERA = pkg_sims.modules['lwMaskERA5'].longitude
latERA = pkg_sims.modules['lwMaskERA5'].latitude

survERA = {}
for key, data in surv.items():
    x, y = surv.parse_key(key)
    if x < 0.0:
        x += 360.0
    xi = getIndex(x, lonERA)
    yi = getIndex(y, latERA)
    survERA.setdefault((xi, yi), []).append(data)

# Average observations within each ERA5 cell
survERA = {
    key: numpy.nanmean(values, axis=0).tolist()
    for key, values in survERA.items()
}

def make_key(lon, lat):
    return f"{lon:.6f},{lat:.6f}"

def parse_key(key):
    return tuple(map(float, key.split(",")))

def getMatrixCERRA(srv, resolution=5500):
    """
    Convert surveillance coordinates into CERRA grid cells.

    Parameters
    ----------
    srv : dict
        Dictionary with keys formatted as "lon,lat".
    crs_cerra : CRS
        Native CERRA projection (e.g. from a CERRA dataset).
    resolution : float
        Grid resolution in metres (default 5500).

    Returns
    -------
    geopandas.GeoDataFrame
        CERRA grid cell polygons in EPSG:4326.
    """

    # Extract coordinates
    coords = numpy.array([parse_key(k) for k in dict.fromkeys(srv)])

    # Create points in geographic coordinates
    gdf = gpd.GeoDataFrame(
        {"key": list(srv.keys()),
         "lon": coords[:, 0],
         "lat": coords[:, 1]},
        geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]),
        crs="EPSG:4326"
    )

    crs_cerra = (
        "+proj=lcc +lat_1=30 +lat_2=65 "
        "+lat_0=50 +lon_0=8 "
        "+R=6371229 +units=m +no_defs"
    )

    # Transform to CERRA projection
    gdf = gdf.to_crs(crs_cerra)

    # Construct square cells centred on each point
    half = resolution / 2

    gdf["geometry"] = [
        box(x - half, y - half, x + half, y + half)
        for x, y in zip(gdf.geometry.x, gdf.geometry.y)
    ]

    # Return geographic polygons
    return gdf.to_crs("EPSG:4326")
