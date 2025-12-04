import json

from ..response import empty_response, returnResponse
from ..functions import isel
import veclim_data_server.pkg_sims as pkg_sims
import veclim_data_server.pkg_surv as pkg_surv

def respond(start_response, kw):
    sims = pkg_sims.modules['Bologna2024']
    surv = pkg_surv.modules['Bologna2024']
    #
    if not 'lon' in kw:
        return returnResponse(start_response, 'Missing argument: lon')
    lon = kw['lon']
    #
    if not 'lat' in kw:
        return returnResponse(start_response, 'Missing argument: lat')
    lat = kw['lat']
    #
    simclm = {
        model: isel(sims.ncs[model].eggs,lon=lon,lat=lat).tolist()
        for model in sims.models
    }
    #
    if not simclm:
        return returnResponse(start_response, empty_response)
    #
    tmp = isel(surv.surv.eggs,lon=lon,lat=lat,ret_df=True)
    survclm = {
        'eggs': tmp.tolist()
    }
    #
    if not survclm:
        return returnResponse(start_response, empty_response)
    #
    ret = {
        'location': {'lon': float(tmp.lon), 'lat': float(tmp.lat)},
        'date': sims.dates.tolist(),
        'models': simclm,
        'surv': survclm
    }
    #
    response_body = json.dumps(ret)
    return returnResponse(start_response, response_body)
