import json

from ..response import empty_response, returnResponse
from ..fun_server import get_decadal
import veclim_data_server.pkg_surv as pkg_surv

def respond(date0,
            date1,
            lon,
            lat,
            timeseries,
            meteo_key,
            sim_key,
            fcast_key,
            risk_key,
            start_response):
    #
    simclm = get_decadal(lon,lat,date0,date1,ts=timeseries)
    if not simclm:
        return returnResponse(start_response, empty_response)
    #
    vec = pkg_surv.modules['albosurv'].presence.search(lon,lat)
    #
    ret = {
        'location': simclm['location'],
        'date': simclm['date'],
        'presence': {
            'albopictus': vec
        }
    }
    #
    if (not simclm['location']['island']) or (not simclm['date']['valid']):
        return returnResponse(start_response, json.dumps(ret))
    #
    ret[meteo_key] = simclm['clm']
    ret[sim_key] = simclm['sim']
    ret[fcast_key] = simclm['fcast']
    ret[risk_key] = simclm['risk']
    #
    if 'surv' in simclm:
        ret['surv-ts'] = simclm['surv']
    #
    response_body = json.dumps(ret)
    return returnResponse(start_response, response_body)
