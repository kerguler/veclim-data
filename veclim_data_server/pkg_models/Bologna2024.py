import json
import numpy
from datetime import date

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
    for model in sims.ensembles:
        simclm[model] = isel(sims.ncs[model].eggs,lon=lon,lat=lat).tolist()
    #
    if not simclm:
        return returnResponse(start_response, empty_response)
    #
    tmp = isel(surv.surv_annual.eggs,lon=lon,lat=lat,ret_df=True)
    survclm = {
        'eggs': tmp.values.tolist()
    }
    #
    if not survclm:
        return returnResponse(start_response, empty_response)
    #
    #dates = {
    #    "date0": f"{date.today().year}-01-01",
    #    "date1": f"{date.today().year}-12-31",
    #    "days": [
    #        0,
    #        51
    #    ],
    #    "valid": 1
    #}
    dates = [f"{date.today().year}-{w:02d}" for w in numpy.arange(1, 53)] 
    #
    ret = {
        'location': {'lon': float(tmp.lon), 'lat': float(tmp.lat)},
        'date': dates,
        'models': simclm,
        'surv': survclm
    }
    #
    response_body = json.dumps(ret)
    return returnResponse(start_response, response_body)
