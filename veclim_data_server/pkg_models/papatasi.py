import json
import numpy

from ..response import empty_response, returnResponse
from ..functions import get_dates, get_clim, getIndex, remove_feb29
import veclim_data_server.pkg_sims as pkg_sims

def get_location(lon, lat, lons, lats):
    papatasi2015 = pkg_sims.modules['papatasi2015']
    #
    if (min(lons) >= 0.0) and (max(lons) >= 180.0) and (lon < 0.0):
        lon += 360.0
    #
    loni = getIndex(lon, lons)
    lati = getIndex(lat, lats)
    #
    lon = lons[loni]
    lat = lats[lati]
    #
    island = papatasi2015.island(loni,lati)
    #
    return {
        'lon': lon,
        'lat': lat,
        'loni': int(loni),
        'lati': int(lati),
        'island': int(island)
    }

def get_papatasi_days(loni, lati, idates, isFeb29):
    papatasi2015 = pkg_sims.modules['papatasi2015']
    #
    return {
        "simL": remove_feb29(papatasi2015.simGERI[:,lati,loni],idates,isFeb29),
        "simH": remove_feb29(papatasi2015.simSTENI[:,lati,loni],idates,isFeb29)
    }

def get_sandfly(lon, lat, date0, date1=False, ts=False):
    papatasi2015 = pkg_sims.modules['papatasi2015']
    #
    dats = get_dates(date0, date1=date1, ts=ts)
    ret = {
        'location': get_location(lon, lat, papatasi2015.longitude, papatasi2015.latitude),
        'date': {key:dats[key] for key in ['date0','date1','days','valid']},
        'clm': {},
        'sim': {},
        'risk': {}
    }
    ret['date']['days'] = ret['date']['days'][[0,-1]].tolist()
    if ((not ret['location']['island']) or 
        (not ret['date']['valid']) or 
        numpy.any([d < 90 for d in dats['days']]) or
        numpy.any([d.year != 2015 for d in dats['dates']])):
        return ret
    #
    ret['sim'] = {
        '2015': get_clim(get_papatasi_days,
                         ret['location']['loni'], 
                         ret['location']['lati'], 
                         dats['days'], 
                         dats['isFeb29'],
                         ts=ts)
    }
    #
    return ret

def respond(start_response, kw):
    if not 'date0' in kw:
        return returnResponse(start_response, 'Missing argument: date0')
    date0 = kw['date0']
    #
    if not 'date1' in kw:
        return returnResponse(start_response, 'Missing argument: date1')
    date1 = kw['date1']
    #
    if not 'lon' in kw:
        return returnResponse(start_response, 'Missing argument: lon')
    lon = kw['lon']
    #
    if not 'lat' in kw:
        return returnResponse(start_response, 'Missing argument: lat')
    lat = kw['lat']
    #
    if not 'timeseries' in kw:
        return returnResponse(start_response, 'Missing argument: timeseries')
    timeseries = kw['timeseries']
    #
    if not 'sim_key' in kw:
        return returnResponse(start_response, 'Missing argument: sim_key')
    sim_key = kw['sim_key']
    #
    simclm = get_sandfly(lon,lat,date0,date1,ts=timeseries)
    if not simclm:
        return returnResponse(start_response, empty_response)
    #
    ret = {
        'location': simclm['location'],
        'date': simclm['date']
    }
    #
    if (not simclm['location']['island']) or (not simclm['date']['valid']):
        return returnResponse(start_response, json.dumps(ret))
    #
    ret[sim_key] = simclm['sim']
    #
    response_body = json.dumps(ret)
    return returnResponse(start_response, response_body)
