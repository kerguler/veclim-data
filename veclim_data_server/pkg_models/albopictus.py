import json
import numpy

from ..response import empty_response, returnResponse
from ..functions import remove_feb29, calc_cut, isel, getIndex, get_dates, get_clim
import veclim_data_server.pkg_surv as pkg_surv
import veclim_data_server.pkg_sims as pkg_sims

def get_location_ERA5(lon, lat):
    lwMaskERA5 = pkg_sims.modules['lwMaskERA5']
    #
    if lon < 0.0:
        lon += 360.0
    #
    loni = getIndex(lon, lwMaskERA5.longitude)
    lati = getIndex(lat, lwMaskERA5.latitude)
    #
    lon = lwMaskERA5.longitude[loni]
    lat = lwMaskERA5.latitude[lati]
    #
    island = lwMaskERA5.island(loni,lati)
    #
    return {
        'lon': lon,
        'lat': lat,
        'loni': int(loni),
        'lati': int(lati),
        'island': int(island)
    }

def get_meteo_days(loni, lati, idates, isFeb29):
    annualERA5 = pkg_sims.modules['annualERA5']
    popdens = pkg_sims.modules['popdens']
    return {
        "photo": remove_feb29(annualERA5.photo[lati,:],idates,isFeb29),
        "atemp": remove_feb29(isel(annualERA5.atemp,lat=lati,lon=loni),idates,isFeb29),
        "atmin": remove_feb29(isel(annualERA5.atmin,lat=lati,lon=loni),idates,isFeb29),
        "atmax": remove_feb29(isel(annualERA5.atmax,lat=lati,lon=loni),idates,isFeb29),
        "rehum": remove_feb29(isel(annualERA5.rehum,lat=lati,lon=loni),idates,isFeb29),
        "precp": remove_feb29(isel(annualERA5.precp,lat=lati,lon=loni),idates,isFeb29),
        "soilw": remove_feb29(isel(annualERA5.soilw,lat=lati,lon=loni),idates,isFeb29),
        "pdens": popdens.pop[lati,loni]
    }

def get_risk(ret):
    labels = [0, 1, 2]
    larva_limit = [-numpy.inf, 0.07, 3.6, numpy.inf]
    adult_limit = [-numpy.inf, 0.7, 3.6, numpy.inf]
    pouts_limit = [-numpy.inf, 1.0, 50.0, numpy.inf]
    iouts_limit = [-numpy.inf, 10.0, 500.0, numpy.inf]
    return {
        "larva": calc_cut(ret["coln2"], larva_limit, labels),
        "adult": calc_cut(ret["colegg"], adult_limit, labels),
        "pouts": calc_cut(ret["pouts"], pouts_limit, labels),
        "iouts": calc_cut(ret["iouts"], iouts_limit, labels)
    }

def get_sim_days(loni, lati, idates, isFeb29):
    annualVectorA = pkg_sims.modules['annualVectorA']
    return {
        "colegg": remove_feb29(isel(annualVectorA.colegg,lat=lati,lon=loni),idates,isFeb29),
        "colK": remove_feb29(isel(annualVectorA.colK,lat=lati,lon=loni),idates,isFeb29),
        "coln2": remove_feb29(isel(annualVectorA.coln2,lat=lati,lon=loni),idates,isFeb29),
        "coln4f": remove_feb29(isel(annualVectorA.coln4f,lat=lati,lon=loni),idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(isel(annualVectorA.pouts,lat=lati,lon=loni),idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(isel(annualVectorA.iouts,lat=lati,lon=loni),idates,isFeb29)]
    }

def get_sim1980_days(loni, lati, idates, isFeb29):
    annualVectorA_1980 = pkg_sims.modules['annualVectorA_1980']
    return {
        "colegg": remove_feb29(isel(annualVectorA_1980.colegg,lat=lati,lon=loni),idates,isFeb29),
        "colK": remove_feb29(isel(annualVectorA_1980.colK,lat=lati,lon=loni),idates,isFeb29),
        "coln2": remove_feb29(isel(annualVectorA_1980.coln2,lat=lati,lon=loni),idates,isFeb29),
        "coln4f": remove_feb29(isel(annualVectorA_1980.coln4f,lat=lati,lon=loni),idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(isel(annualVectorA_1980.pouts,lat=lati,lon=loni),idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(isel(annualVectorA_1980.iouts,lat=lati,lon=loni),idates,isFeb29)]
    }

def get_nasa_ssp245_days(loni, lati, idates, isFeb29):
    annualNASA = pkg_sims.modules['annualNASA']
    ssp = 'ssp245'
    return {
        "colegg": remove_feb29(annualNASA.colegg[ssp][lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualNASA.pouts[ssp][lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualNASA.iouts[ssp][lati,loni,:],idates,isFeb29)]
    }

def get_nasa_ssp585_days(loni, lati, idates, isFeb29):
    annualNASA = pkg_sims.modules['annualNASA']
    ssp = 'ssp585'
    return {
        "colegg": remove_feb29(annualNASA.colegg[ssp][lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualNASA.pouts[ssp][lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualNASA.iouts[ssp][lati,loni,:],idates,isFeb29)]
    }

def get_fcast_days(loni, lati, date0, date1):
    forecastECMWF = pkg_sims.modules['forecastECMWF']
    xr = (forecastECMWF.dates >= numpy.datetime64(date0)) & (forecastECMWF.dates <= numpy.datetime64(date1))
    if not any(xr):
        return {}
    return {
        "colegg": numpy.nan_to_num(isel(forecastECMWF.colegg,lat=lati,lon=loni)[xr],nan=0.0),
        "colK": numpy.nan_to_num(isel(forecastECMWF.colK,lat=lati,lon=loni)[xr],nan=0.0),
        "coln2": numpy.nan_to_num(isel(forecastECMWF.coln2,lat=lati,lon=loni)[xr],nan=0.0),
        "coln4f": numpy.nan_to_num(isel(forecastECMWF.coln4f,lat=lati,lon=loni)[xr],nan=0.0),
        "pouts": numpy.nan_to_num(isel(forecastECMWF.pouts,lat=lati,lon=loni)[xr[:-60]],nan=0.0)*100.0,
        "iouts": numpy.nan_to_num(isel(forecastECMWF.iouts,lat=lati,lon=loni)[xr[:-60]],nan=0.0)*4000.0
    }

def get_surv(lon, lat, idates, isFeb29):
    albosurv = pkg_surv.modules['albosurv']
    return {
        key: [] if len(value)==0 else remove_feb29(value,idates,isFeb29,fill=None)
        for key,value in albosurv.getSurv(lon,lat).items()
    }

def get_decadal(lon, lat, date0, date1=False, ts=False):
    annualVectorA = pkg_sims.modules['annualVectorA']
    forecastECMWF = pkg_sims.modules['forecastECMWF']
    annualNASA = pkg_sims.modules['annualNASA']
    #
    dats = get_dates(date0, date1=date1, ts=ts)
    ret = {
        'location': get_location_ERA5(lon, lat),
        'date': {key:dats[key] for key in ['date0','date1','valid']},
        'clm': {},
        'sim': {},
        'risk': {}
    }
    if (not ret['location']['island']) or (not ret['date']['valid']):
        return ret
    #
    ret['clm'] = {
        '2010-2019': get_clim(get_meteo_days,
                              ret['location']['loni'], 
                              ret['location']['lati'], 
                              dats['days'], 
                              dats['isFeb29'],
                              ts=ts)
    }
    #
    acc = annualVectorA.accuracy(ret['location']['lati'],ret['location']['loni'])
    ret['sim'] = {
        '1980-1989': get_clim(get_sim1980_days, 
                              ret['location']['loni'], 
                              ret['location']['lati'], 
                              dats['days'], 
                              dats['isFeb29'],
                              ts=ts),
        '2010-2019': get_clim(get_sim_days, 
                              ret['location']['loni'], 
                              ret['location']['lati'], 
                              dats['days'], 
                              dats['isFeb29'],
                              ts=ts)
    }
    ret['sim']['2010-2019']['accuracy'] = acc
    #
    dts = forecastECMWF.getOverlap()
    ret['fcast'] = {
        'ecmwf': get_clim(get_fcast_days, 
                              ret['location']['loni'], 
                              ret['location']['lati'], 
                              date0, 
                              date1,
                              ts=ts),
        'nasa': {
            'ssp245': get_clim(get_nasa_ssp245_days, 
                              ret['location']['loni'], 
                              ret['location']['lati'], 
                              dats['days'], 
                              dats['isFeb29'],
                              ts=ts),
            'ssp585': get_clim(get_nasa_ssp585_days, 
                              ret['location']['loni'], 
                              ret['location']['lati'], 
                              dats['days'], 
                              dats['isFeb29'],
                              ts=ts)
        }
    }
    ret['fcast']['ecmwf']['overlap'] = dts
    ret['fcast']['nasa']['decade'] = annualNASA.decade
    ret['fcast']['nasa']['models'] = annualNASA.models
    #
    tmp = get_risk(ret['sim']['2010-2019'])
    ret['risk'] = {
        '2010-2019': {
            'adult': tmp['adult'],
            'larva': tmp['larva'],
            'pouts': tmp['pouts'],
            'iouts': tmp['iouts'],
            'accuracy': acc
        }
    }
    #
    if ts:
        ret['surv'] = get_surv(ret['location']['lon'],
                               ret['location']['lat'],
                               dats['days'],
                               dats['isFeb29'])
    #
    return ret

def respond(start_response, kw):
    #
    print("LOG: Respondidng to",kw,flush=True)
    #
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
    if not 'meteo_key' in kw:
        return returnResponse(start_response, 'Missing argument: meteo_key')
    meteo_key = kw['meteo_key']
    #
    if not 'sim_key' in kw:
        return returnResponse(start_response, 'Missing argument: sim_key')
    sim_key = kw['sim_key']
    #
    if not 'fcast_key' in kw:
        return returnResponse(start_response, 'Missing argument: fcast_key')
    fcast_key = kw['fcast_key']
    #
    if not 'risk_key' in kw:
        return returnResponse(start_response, 'Missing argument: risk_key')
    risk_key = kw['risk_key']
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
