import json
import numpy

from ..response import returnResponse
from ..functions import get_dates, remove_feb29
import veclim_data_server.pkg_sims as pkg_sims
import veclim_data_server.pkg_surv as pkg_surv

def get_sandfly(lon, lat, date0, date1=False, ts=False):
    papatasi = pkg_sims.modules['papatasi_V2511A_PRT']
    surv = pkg_surv.modules['papatasi_V2511A_PRT']
    ret = {
        'location': {
            'lon': lon,
            'lat': lat,
            'pid': None,
            'name': None,
            'island': 0
        },
        'date': {}
    }
    #
    prop = papatasi.getPolyProp(papatasi.db.polys,lon,lat)
    if prop is None:
        return ret
    pid = prop['Official_Co']
    name = prop['Official_Na']
    ret['location']['pid'] = pid
    ret['location']['name'] = name
    ret['location']['island'] = 1
    #
    dats = get_dates(date0, date1=date1, ts=ts)
    ret['date'] = {key:dats[key] for key in ['date0','date1','days','valid']}
    ret['date']['days'] = ret['date']['days'][[0,-1]].tolist()
    if not ret['date']['valid']:
        return ret
    #
    val = remove_feb29(papatasi.db.mat[papatasi.db.var].sel(poly=pid),
                       dats['days'], dats['isFeb29'], fill=None, mean=True)
    #
    # times: 0..364
    up = papatasi.up_times[pid]
    down = papatasi.down_times[pid]
    peak_up = papatasi.peak_up_times[pid]
    peak_down = papatasi.peak_down_times[pid]
    #
    risk = papatasi.classify_days(dats['days']+1, 
                                  up, 
                                  down, 
                                  peak_up, 
                                  peak_down)
    #
    srv = []
    if pid in surv.surv:
        srv = [surv.surv[pid][doy] for doy in dats['days']]
    #
    ret['sim'] = {
        'V2511A_PRT': {
            papatasi.db.var: val
        }
    }
    #
    ret['risk'] = {
        'V2511A_PRT': {
            'up': up,
            'down': down,
            'peak_up': peak_up,
            'peak_down': peak_down,
            'risk': risk
        }
    }
    #
    ret['surv'] = {
        'adult_norm': srv
    }
    #
    return ret

def respond(start_response, kw):
    if not 'date0' in kw:
        return returnResponse(start_response, 'Missing argument: date0')
    date0 = kw['date0']
    #
    if date0 is None:
        return returnResponse(start_response, 'Missing argument: date0')
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
    if not 'risk_key' in kw:
        return returnResponse(start_response, 'Missing argument: risk_key')
    risk_key = kw['risk_key']
    #
    simclm = get_sandfly(lon,lat,date0,date1,ts=timeseries)
    if ((not simclm) or
        (not simclm['location']['island']) or 
        (not simclm['date']['valid'])):
        return returnResponse(start_response, json.dumps(simclm))
    #
    ret = {
        'location': simclm['location'],
        'date': simclm['date']
    }
    ret[sim_key] = simclm['sim']
    ret[risk_key] = simclm['risk']
    #
    if 'surv' in simclm:
        ret['surv-ts'] = simclm['surv']
    #
    response_body = json.dumps(ret)
    return returnResponse(start_response, response_body)
