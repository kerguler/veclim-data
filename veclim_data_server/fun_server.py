import numpy
import pandas

import veclim_data_server.pkg_sims as pkg_sims
import veclim_data_server.pkg_tiles as pkg_tiles
import veclim_data_server.pkg_surv as pkg_surv

from veclim_data_server.functions import is_leap_year, getIndex

from datetime import datetime, timedelta
Feb29 = datetime(2020,2,29).timetuple().tm_yday

veclist = ['albopictus', 'papatasi']

not_clim_keys = [
    "inv",
    "lon",
    "lat",
    "day0",
    "day1",
]

def remove_feb29(vec,days,isFeb29,fill=0.0):
    tmp = vec[days-1]
    if len(isFeb29) == 0:
        if fill == None:
            return [None if numpy.isnan(a) else a for a in tmp]
        else:
            return numpy.nan_to_num(tmp,nan=fill).tolist()
    tmp[isFeb29+1] = 0.5*(tmp[isFeb29]+tmp[isFeb29+1])
    if fill == None:
        return [None if numpy.isnan(a) else a for a in numpy.delete(tmp,isFeb29)]
    else:
        return numpy.nan_to_num(numpy.delete(tmp,isFeb29),nan=fill).tolist()

def calc_list(clms):
    return {
        key: numpy.hstack([clm[key] for clm in clms]).tolist()
        for key in clms[0] if not key in not_clim_keys
    }

def calc_mean(clms):
    return {
        key: float(numpy.nanmean(numpy.hstack([clm[key] for clm in clms])))
        for key in clms[0] if not key in not_clim_keys
    }

def get_clim(get_days, loni, lati, pr0, pr1, ts=False):
    if ts:
        calc_fun = calc_list
    else:
        calc_fun = calc_mean
    #
    ret = calc_fun([get_days(loni,lati,pr0,pr1)])
    #
    return ret

def get_dates(date0, date1=False, ts=False):
    valid = True
    if (not date1) or (date1 < date0):
        date1 = date0
        valid = False
    #
    isFeb29 = []
    idates = []
    ddates = []
    while date0 <= date1:
        ddates.append(date0)
        idates.append(date0.timetuple().tm_yday-1 if is_leap_year(date0.year) and date0.timetuple().tm_yday>Feb29 else date0.timetuple().tm_yday)
        if date0.month==2 and date0.day==29:
            isFeb29.append(len(ddates)-1)
            ddates.append(date0)
            idates.append(Feb29-1)
        date0 += timedelta(days=1)
    #
    ddates = numpy.array(ddates)
    idates = numpy.array(idates)
    isFeb29 = numpy.array(isFeb29)
    #
    return {
        "dates": ddates,
        "days": idates,
        "isFeb29": isFeb29,
        "date0": ddates[0].strftime("%Y-%m-%d"),
        "date1": ddates[-1].strftime("%Y-%m-%d"),
        "valid": int(valid)
    }

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
        "atemp": remove_feb29(annualERA5.atemp[lati,loni,:],idates,isFeb29),
        "atmin": remove_feb29(annualERA5.atmin[lati,loni,:],idates,isFeb29),
        "atmax": remove_feb29(annualERA5.atmax[lati,loni,:],idates,isFeb29),
        "rehum": remove_feb29(annualERA5.rehum[lati,loni,:],idates,isFeb29),
        "precp": remove_feb29(annualERA5.precp[lati,loni,:],idates,isFeb29),
        "soilw": remove_feb29(annualERA5.soilw[lati,loni,:],idates,isFeb29),
        "pdens": popdens.pop[lati,loni]
    }

def calc_cut(vec,lim,lab):
    if hasattr(vec, '__iter__'):
        return pandas.cut(vec, bins=lim, include_lowest=False, right=False, labels=lab).tolist()
    return lab[numpy.where(numpy.array(lim) > vec)[0][0] - 1]

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
        "colegg": remove_feb29(annualVectorA.colegg[lati,loni,:],idates,isFeb29),
        "colK": remove_feb29(annualVectorA.colK[lati,loni,:],idates,isFeb29),
        "coln2": remove_feb29(annualVectorA.coln2[lati,loni,:],idates,isFeb29),
        "coln4f": remove_feb29(annualVectorA.coln4f[lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualVectorA.pouts[lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualVectorA.iouts[lati,loni,:],idates,isFeb29)]
    }

def get_sim1980_days(loni, lati, idates, isFeb29):
    annualVectorA_1980 = pkg_sims.modules['annualVectorA_1980']
    return {
        "colegg": remove_feb29(annualVectorA_1980.colegg[lati,loni,:],idates,isFeb29),
        "colK": remove_feb29(annualVectorA_1980.colK[lati,loni,:],idates,isFeb29),
        "coln2": remove_feb29(annualVectorA_1980.coln2[lati,loni,:],idates,isFeb29),
        "coln4f": remove_feb29(annualVectorA_1980.coln4f[lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualVectorA_1980.pouts[lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualVectorA_1980.iouts[lati,loni,:],idates,isFeb29)]
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
        "colegg": numpy.nan_to_num(forecastECMWF.colegg[lati,loni,xr],nan=0.0),
        "colK": numpy.nan_to_num(forecastECMWF.colK[lati,loni,xr],nan=0.0),
        "coln2": numpy.nan_to_num(forecastECMWF.coln2[lati,loni,xr],nan=0.0),
        "coln4f": numpy.nan_to_num(forecastECMWF.coln4f[lati,loni,xr],nan=0.0),
        "pouts": numpy.nan_to_num(forecastECMWF.pouts[lati,loni,xr[:-60]],nan=0.0)*100.0,
        "iouts": numpy.nan_to_num(forecastECMWF.iouts[lati,loni,xr[:-60]],nan=0.0)*4000.0
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