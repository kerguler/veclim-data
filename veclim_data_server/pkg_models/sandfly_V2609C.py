import json
import numpy
import pandas

from ..response import empty_response, returnResponse
from ..functions import get_dates, get_clim, getIndex, isel_feb29, remove_feb29
import veclim_data_server.pkg_sims as pkg_sims
import veclim_data_server.pkg_surv as pkg_surv

def get_surv_SHP(pid, dates):
    surv = pkg_surv.modules['sandfly_EDENext'].survPRT
    #
    srv = []
    if pid in surv:
        srv = [surv[pid][doy] for doy in dates['days']]
    #
    return {
        'adult_norm': srv
    }

def get_surv_ERA(loni, lati, dates):
    surv = pkg_surv.modules['sandfly_EDENext'].survERA
    #
    data = surv.get((loni, lati))
    #
    if data is None:
        return [numpy.nan] * len(dates['days'])
    #
    srv = [data[doy] for doy in dates['days']]
    #
    return {
        'adult_norm': srv
    }

def get_surv_CERRA(lon, lat, dates):
    surv = pkg_surv.modules['sandfly_EDENext']
    #
    data = surv.surv.get(surv.make_key(lon, lat))
    #
    if data is None:
        return [numpy.nan] * len(dates['days'])
    #
    srv = [data[doy] for doy in dates['days']]
    #
    return {
        'adult_norm': srv
    }

def get_location_SHP(lon, lat):
    papatasi = pkg_sims.modules['papatasi_V2609C_PRT']
    #
    loc = {
        'lon': lon,
        'lat': lat,
        'pid': None,
        'name': None,
        'island': 0
    }
    #
    prop = papatasi.getPolyProp(papatasi.db.polys,lon,lat)
    if prop is None:
        return loc
    pid = prop['Official_Co']
    name = prop['Official_Na']
    loc['pid'] = pid
    loc['name'] = name
    loc['island'] = 1
    #
    return loc

def get_location_CERRA(lon, lat):
    grid = pkg_sims.modules['gridCERRA'].grid
    #
    grd = grid.geti(lon,lat)
    #
    island = grid.ds.lsm[grd['y'],grd['x']]
    #
    return {
        'lon': grd['lon'],
        'lat': grd['lat'],
        'loni': grd['x'],
        'lati': grd['y'],
        'island': int(island > 0.1)
    }

def get_location_ERA(lon, lat):
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

def get_sandfly_ERA(loni, lati, idates, isFeb29, ts=True):
    papatasi = pkg_sims.modules['papatasi_V2609C']
    #
    return {
        "newegg": isel_feb29(papatasi.newegg['newegg'],loni=loni,lati=lati,idates=idates,isFeb29=isFeb29),
        "newegg_ssp585": isel_feb29(papatasi.newegg_ssp585['newegg'],loni=loni,lati=lati,idates=idates,isFeb29=isFeb29)
    }

def get_sandfly_CERRA(loni, lati, idates, isFeb29, ts=True):
    papatasi = pkg_sims.modules['papatasi_V2609C']
    #
    return {
        "newegg": isel_feb29(papatasi.newegg_CERRA['newegg'],x=loni,y=lati,idates=idates,isFeb29=isFeb29),
    }

def get_sandfly_SHP(pid, dates):
    papatasi = pkg_sims.modules['papatasi_V2609C_PRT']
    #
    val = remove_feb29(papatasi.db.mat[papatasi.db.var].sel(poly=pid), dates['days'], dates['isFeb29'], fill=None, mean=True)
    #
    sim = {
        'V2609C_PRT': {
            papatasi.db.var: val
        }
    }
    #
    return sim

def get_sandfly_risk_SHP(pid, dates):
    papatasi = pkg_sims.modules['papatasi_V2609C_PRT']
    #
    ret = {
    }
    #
    # times: 0..364
    up = papatasi.up_times[pid]
    down = papatasi.down_times[pid]
    peak_up = papatasi.peak_up_times[pid]
    peak_down = papatasi.peak_down_times[pid]
    #
    risk = papatasi.classify_days(pandas.to_datetime(dates['dates']), 
                                  up, 
                                  down, 
                                  peak_up, 
                                  peak_down)
    #
    rsk = {
        'V2609C_PRT': {
            'up': [t.strftime('%Y-%m-%d') for t in up],
            'down': [t.strftime('%Y-%m-%d') for t in down],
            'peak_up': [t.strftime('%Y-%m-%d') for t in peak_up],
            'peak_down': [t.strftime('%Y-%m-%d') for t in peak_down],
            'risk': risk
        }
    }
    #
    return rsk

def get_sandfly_alert_SHP(pid):
    papatasi = pkg_sims.modules['papatasi_V2609C_PRT']
    #
    notify = papatasi.getWarning(pid)
    #
    return notify

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
    if not 'risk_key' in kw:
        return returnResponse(start_response, 'Missing argument: risk_key')
    risk_key = kw['risk_key']
    #
    rawdates = get_dates(date0, date1=date1, ts=timeseries)
    dates =  {key:rawdates[key] for key in ['date0','date1','days','valid']}
    dates['days'] = dates['days'][[0,-1]].tolist()
    if not dates['valid']:
        return returnResponse(start_response, empty_response)
    #
    locSHP = get_location_SHP(lon, lat)
    locERA = get_location_ERA(lon, lat)
    locCERRA = get_location_CERRA(lon, lat)
    #
    if (not locERA['island'] and 
        not locCERRA['island']):
        return returnResponse(start_response, json.dumps(ret))
    #
    simSHP = {}
    riskSHP = {}
    alertSHP = {}
    survSHP = {}
    if not locSHP['pid'] is None:
        simSHP = get_sandfly_SHP(locSHP['pid'], 
                                 rawdates)
        riskSHP = get_sandfly_risk_SHP(locSHP['pid'], 
                                       rawdates)
        alertSHP = get_sandfly_alert_SHP(locSHP['pid'])
        survSHP = get_surv_SHP(locSHP['pid'], 
                               rawdates)
    #
    simERA = {
        'V2609C': get_clim(get_sandfly_ERA,
                           locERA['loni'], 
                           locERA['lati'], 
                           rawdates['days'], 
                           rawdates['isFeb29'],
                           ts=timeseries)
    }
    survERA = get_surv_ERA(locERA['loni'], 
                           locERA['lati'], 
                           rawdates)
    #
    simCERRA = {
        'V2609C': get_clim(get_sandfly_CERRA,
                           locCERRA['loni'], 
                           locCERRA['lati'], 
                           rawdates['days'], 
                           rawdates['isFeb29'],
                           ts=timeseries)
    }
    survCERRA = get_surv_CERRA(locCERRA['lon'], 
                               locCERRA['lat'], 
                               rawdates)
    #
    ret = {
        'date': dates,
        'SHP': {
            'location': locSHP,
        },
        'ERA': {
            'location': locERA,
        },
        'CERRA': {
            'location': locCERRA,
        }
    }
    #
    if simSHP:
        ret['SHP'][sim_key] = simSHP
    if survSHP:
        ret['SHP']['surv-ts'] = survSHP
    if riskSHP:
        ret['SHP'][risk_key] = riskSHP
    if alertSHP:
        ret['SHP']['alert-ts'] = alertSHP
    #
    if simERA:
        ret['ERA'][sim_key] = simERA
    if survERA:
        ret['ERA']['surv-ts'] = survERA
    #
    if simCERRA:
        ret['CERRA'][sim_key] = simCERRA
    if survCERRA:
        ret['CERRA']['surv-ts'] = survCERRA
    #
    response_body = json.dumps(ret)
    return returnResponse(start_response, response_body)
