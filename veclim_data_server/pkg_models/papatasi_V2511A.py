import json
import numpy

from ..response import empty_response, returnResponse
from ..functions import get_dates, get_clim, isel, remove_feb29
import veclim_data_server.pkg_sims as pkg_sims

def get_location(lon, lat):
    papatasi = pkg_sims.modules['papatasi_V2511A']
    #
    grd = papatasi.grid.geti(lon,lat)
    #
    island = papatasi.grid.ds.lsm[grd['y'],grd['x']]
    #
    return {
        'lon': grd['lon'],
        'lat': grd['lat'],
        'loni': grd['x'],
        'lati': grd['y'],
        'island': int(island > 0.1)
    }

def get_papatasi_days(loni, lati, idates, isFeb29, ts=True):
    papatasi = pkg_sims.modules['papatasi_V2511A']
    #
    return {
        "simL": remove_feb29(isel(papatasi.female_lo['female_lo'],x=loni,y=lati),idates,isFeb29),
        "simH": remove_feb29(isel(papatasi.female_hi['female_hi'],x=loni,y=lati),idates,isFeb29),
        "simM": remove_feb29(isel(papatasi.female_mn['female_mn'],x=loni,y=lati),idates,isFeb29)
    }

def get_sandfly(lon, lat, date0, date1=False, ts=False):
    dats = get_dates(date0, date1=date1, ts=ts)
    ret = {
        'location': get_location(lon, lat),
        'date': {key:dats[key] for key in ['date0','date1','days','valid']},
        'clm': {},
        'sim': {},
        'risk': {}
    }
    ret['date']['days'] = ret['date']['days'][[0,-1]].tolist()
    if ((not ret['location']['island']) or 
        (not ret['date']['valid'])):
        return ret
    #
    ret['sim'] = {
        'V2511A': get_clim(get_papatasi_days,
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
