import numpy
import pandas

import fun_clim
import fun_surv
import fun_colors
import fun_tiles

Feb29 = fun_clim.datetime(2020,2,29).timetuple().tm_yday

veclist = ['albopictus', 'papatasi']

not_clim_keys = [
    "inv",
    "lon",
    "lat",
    "day0",
    "day1",
]

lwMaskERA5 = False
annualERA5 = False
annualVectorA = False
annualVectorA_1980 = False
annualNASA = False
forecastECMWF = False
papatasi2015 = False
popdens = False
presence = False
vabun = False
albosurv = False

def load_forecast_var(reload=False):
    global forecastECMWF
    if reload or (not forecastECMWF):
        print("Loading ECMWF forecast...",flush=True)
        forecastECMWF = fun_clim.forecastECMWF()

def load_global_var():
    global lwMaskERA5
    global annualERA5
    global annualVectorA
    global annualVectorA_1980
    global annualNASA
    global papatasi2015
    global popdens
    global presence
    global vabun
    global albosurv
    #
    if not lwMaskERA5:
        print("Loading ERA5 land/water mask...",flush=True)
        lwMaskERA5 = fun_clim.lwMaskERA5()
    #
    if not annualERA5:
        print("Loading annual ERA5...",flush=True)
        annualERA5 = fun_clim.annualERA5()
    #
    if not annualVectorA:
        print("Loading VectorA...",flush=True)
        annualVectorA = fun_clim.annualVectorA()
    #
    if not annualVectorA_1980:
        print("Loading VectorA_1980...",flush=True)
        annualVectorA_1980 = fun_clim.annualVectorA_1980()
    #
    if not annualNASA:
        print("Loading NASA...",flush=True)
        annualNASA = fun_clim.annualNASA()
    #
    if not papatasi2015:
        print("Loading papatasi2015...",flush=True)
        papatasi2015 = fun_clim.papatasi2015()
    #
    if not popdens:
        print("Loading popdens...",flush=True)
        popdens = fun_clim.popDens()
    #
    if not presence:
        print("Loading presence...",flush=True)
        presence = fun_surv.presenceAlbopictus()
    #
    if not vabun:
        print("Loading VectAbundance...",flush=True)
        vabun = fun_surv.VectAbundance()
    #
    if not albosurv:
        print("Loading multiple presence datasets...",flush=True)
        albosurv = fun_surv.albosurv()
    #
    load_forecast_var(reload=False)

def remove_feb29(vec,days,isFeb29):
    tmp = vec[days-1]
    if len(isFeb29) == 0:
        return tmp
    tmp[isFeb29+1] = 0.5*(tmp[isFeb29]+tmp[isFeb29+1])
    return numpy.delete(tmp,isFeb29).tolist()

def remove2_feb29(mat,days,isFeb29):
    tmp = mat[:,days-1]
    if len(isFeb29) == 0:
        return tmp
    tmp[:,isFeb29+1] = 0.5*(tmp[:,isFeb29]+tmp[:,isFeb29+1])
    return numpy.delete(tmp,isFeb29,axis=1).tolist()

def remove3_feb29(mat,days,isFeb29,tolist=True):
    tmp = mat[:,:,days-1]
    if len(isFeb29) == 0:
        return tmp
    tmp[:,:,isFeb29+1] = 0.5*(tmp[:,:,isFeb29]+tmp[:,:,isFeb29+1])
    tmp = numpy.delete(tmp,isFeb29,axis=2)
    return tmp.tolist() if tolist else tmp


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
        idates.append(date0.timetuple().tm_yday-1 if fun_clim.is_leap_year(date0.year) and date0.timetuple().tm_yday>Feb29 else date0.timetuple().tm_yday)
        if date0.month==2 and date0.day==29:
            isFeb29.append(len(ddates)-1)
            ddates.append(date0)
            idates.append(Feb29-1)
        date0 += fun_clim.timedelta(days=1)
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


tile_dat = {}

def load_tiles_dates(v_label,date0=None,date1=None):
    global tile_dat
    #
    load_global_var()
    #
    if ((date0 == None) or (date1 == None)):
        return {}
    #
    dt = get_dates(date0, date1=date1, ts=False)
    #
    if v_label == 'colegg_dates':
        cmap = fun_colors.cmaps['FuzzyLocV6']
        tile_dat['colegg_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](remove3_feb29(annualVectorA.colegg,dt['days'],dt['isFeb29'],tolist=False)[:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    elif v_label == 'larva_dates':
        cmap = fun_colors.cmaps['larva']
        tile_dat['larva_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](remove3_feb29(annualVectorA.coln2,dt['days'],dt['isFeb29'],tolist=False)[:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    elif v_label == 'chikv_iouts_dates':
        cmap = fun_colors.cmaps['iouts']
        tile_dat['chikv_iouts_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](remove3_feb29(annualVectorA.iouts,dt['days'],dt['isFeb29'],tolist=False)[:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    elif v_label == 'chikv_pouts_dates':
        cmap = fun_colors.cmaps['pouts']
        tile_dat['chikv_pouts_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](remove3_feb29(annualVectorA.pouts,dt['days'],dt['isFeb29'],tolist=False)[:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    #
    date0 = pandas.to_datetime(date0)
    date1 = pandas.to_datetime(date1)
    #
    if v_label == 'colegg_fcast_dates':
        dt = (forecastECMWF.dates >= date0) & (forecastECMWF.dates <= date1)
        if numpy.abs(numpy.sum(dt) - (date1-date0).days) > 7:
            return {'error': "Forecast dates do not match the request!"}
        #
        cmap = fun_colors.cmaps['FuzzyLocV6']
        tile_dat['colegg_fcast_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](forecastECMWF.colegg[:-1,:,dt]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    elif v_label == 'larva_fcast_dates':
        dt = (forecastECMWF.dates >= date0) & (forecastECMWF.dates <= date1)
        if numpy.abs(numpy.sum(dt) - (date1-date0).days) > 7:
            return {'error': "Forecast dates do not match the request!"}
        #
        cmap = fun_colors.cmaps['larva']
        tile_dat['larva_fcast_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](forecastECMWF.coln2[:-1,:,dt]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    elif v_label == 'chikv_iouts_fcast_dates':
        dt = (forecastECMWF.idates >= date0) & (forecastECMWF.idates <= date1)
        if numpy.abs(numpy.sum(dt) - (date1-date0).days) > 7:
            return {'error': "Forecast dates do not match the request!"}
        #
        cmap = fun_colors.cmaps['iouts']
        tile_dat['chikv_iouts_fcast_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](forecastECMWF.iouts[:-1,:,dt]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    elif v_label == 'chikv_pouts_fcast_dates':
        dt = (forecastECMWF.idates >= date0) & (forecastECMWF.idates <= date1)
        if numpy.abs(numpy.sum(dt) - (date1-date0).days) > 7:
            return {'error': "Forecast dates do not match the request!"}
        #
        cmap = fun_colors.cmaps['pouts']
        tile_dat['chikv_pouts_fcast_dates'] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](forecastECMWF.pouts[:-1,:,dt]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': '',
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    # 
    return {}

def load_tiles():
    global tile_dat
    #
    load_global_var()
    #
    if 'colegg' not in tile_dat:
        print("Loading tiles: colegg...",flush=True)
        cmap = fun_colors.cmaps['FuzzyLocV6']
        tile_dat['colegg'] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA.colegg[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'colegg',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if 'larva' not in tile_dat:
        print("Loading tiles: larva...",flush=True)
        cmap = fun_colors.cmaps['larva']
        tile_dat['larva'] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA.coln2[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'larva',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if 'chikv_iouts' not in tile_dat:
        print("Loading tiles: chikv_iouts...",flush=True)
        cmap = fun_colors.cmaps['iouts']
        tile_dat['chikv_iouts'] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA.iouts[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'chikv_iouts',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if 'chikv_pouts' not in tile_dat:
        print("Loading tiles: chikv_pouts...",flush=True)
        cmap = fun_colors.cmaps['pouts']
        tile_dat['chikv_pouts'] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA.pouts[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'chikv_pouts',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if "colegg_1980" not in tile_dat:
        print("Loading tiles: colegg_1980...",flush=True)
        cmap = fun_colors.cmaps['FuzzyLocV6']
        tile_dat["colegg_1980"] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA_1980.colegg[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': "colegg_1980",
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    if "chikv_iouts_1980" not in tile_dat:
        print("Loading tiles: chikv_iouts_1980...",flush=True)
        cmap = fun_colors.cmaps['iouts']
        tile_dat["chikv_iouts_1980"] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA_1980.iouts[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': "chikv_iouts_1980",
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    if "chikv_pouts_1980" not in tile_dat:
        print("Loading tiles: chikv_pouts_1980...",flush=True)
        cmap = fun_colors.cmaps['pouts']
        tile_dat["chikv_pouts_1980"] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA_1980.pouts[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': "chikv_pouts_1980",
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if "diff_colegg_1980" not in tile_dat:
        print("Loading tiles: diff_colegg_1980...",flush=True)
        cmap = fun_colors.cmaps['diff_colegg']
        # Based on FuzzyLocV6
        tmp = numpy.log2(numpy.nanmean(annualVectorA.colegg[:-1,:,:],axis=2)/5.0)
        tmp[tmp<-4] = -4
        tmp[tmp>4] = 4
        tmpf = numpy.log2(numpy.nanmean(annualVectorA_1980.colegg[:-1,:,:],axis=2)/5.0)
        tmpf[tmpf<-4] = -4
        tmpf[tmpf>4] = 4
        tmp = tmp-tmpf
        # 
        tile_dat["diff_colegg_1980"] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](tmp),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': "diff_colegg_1980",
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    if "diff_iouts_1980" not in tile_dat:
        print("Loading tiles: diff_iouts_1980...",flush=True)
        cmap = fun_colors.cmaps['diff_iouts']
        tile_dat["diff_iouts_1980"] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA.iouts[:-1,:,:]-annualVectorA_1980.iouts[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': "diff_iouts_1980",
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    if "diff_pouts_1980" not in tile_dat:
        print("Loading tiles: diff_pouts_1980...",flush=True)
        cmap = fun_colors.cmaps['diff_pouts']
        tile_dat["diff_pouts_1980"] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](annualVectorA.pouts[:-1,:,:]-annualVectorA_1980.pouts[:-1,:,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': "diff_pouts_1980",
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    for ssp in ['ssp245','ssp585']:
        if "colegg_%s" %ssp not in tile_dat:
            print("Loading tiles: colegg_%s..." %ssp,flush=True)
            cmap = fun_colors.cmaps['FuzzyLocV6']
            tile_dat["colegg_%s" %ssp] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](annualNASA.colegg[ssp][:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': "colegg_%s" %ssp,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
        if "iouts_%s" %ssp not in tile_dat:
            print("Loading tiles: iouts_%s..." %ssp,flush=True)
            cmap = fun_colors.cmaps['iouts']
            tile_dat["iouts_%s" %ssp] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](annualNASA.iouts[ssp][:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': "iouts_%s" %ssp,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
        if "pouts_%s" %ssp not in tile_dat:
            print("Loading tiles: pouts_%s..." %ssp,flush=True)
            cmap = fun_colors.cmaps['pouts']
            tile_dat["pouts_%s" %ssp] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](annualNASA.pouts[ssp][:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': "pouts_%s" %ssp,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
        if "diff_colegg_%s" %ssp not in tile_dat:
            print("Loading tiles: diff_colegg_%s..." %ssp,flush=True)
            cmap = fun_colors.cmaps['diff_colegg']
            # Based on FuzzyLocV6
            tmp = numpy.log2(numpy.nanmean(annualNASA.colegg[ssp][:-1,:,:],axis=2)/5.0)
            tmp[tmp<-4] = -4
            tmp[tmp>4] = 4
            tmpf = numpy.log2(numpy.nanmean(annualVectorA.colegg[:-1,:,:],axis=2)/5.0)
            tmpf[tmpf<-4] = -4
            tmpf[tmpf>4] = 4
            tmp = tmp-tmpf
            # 
            tile_dat["diff_colegg_%s" %ssp] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](tmp),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': "diff_colegg_%s" %ssp,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
        if "diff_iouts_%s" %ssp not in tile_dat:
            print("Loading tiles: diff_iouts_%s..." %ssp,flush=True)
            cmap = fun_colors.cmaps['diff_iouts']
            tile_dat["diff_iouts_%s" %ssp] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](annualNASA.iouts[ssp][:-1,:,:]-annualVectorA.iouts[:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': "diff_iouts_%s" %ssp,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
        if "diff_pouts_%s" %ssp not in tile_dat:
            print("Loading tiles: diff_pouts_%s..." %ssp,flush=True)
            cmap = fun_colors.cmaps['diff_pouts']
            tile_dat["diff_pouts_%s" %ssp] = {
                'fun': fun_tiles.getTiles,
                'dat': cmap['tran'](annualNASA.pouts[ssp][:-1,:,:]-annualVectorA.pouts[:-1,:,:]),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': "diff_pouts_%s" %ssp,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }
    #
    if 'pop' not in tile_dat:
        print("Loading tiles: pop...",flush=True)
        cmap = fun_colors.cmaps['viridis.pop']
        tile_dat['pop'] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](popdens.pop[:-1,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'pop_viridis_pop',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if 'presence' not in tile_dat:
        print("Loading tiles: presence...",flush=True)
        cmap = fun_colors.cmaps['presence']
        tile_dat['presence'] = {
            'fun': fun_tiles.getTiles,
            'dat': cmap['tran'](presence.matrix[:-1,:]),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'presence',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if 'vabun_v015' not in tile_dat:
        print("Loading tiles: vabun_v015...",flush=True)
        cmap = fun_colors.cmaps['vabun']
        tile_dat['vabun_v015'] = {
            'fun': fun_tiles.getShpTiles,
            'dat': cmap['tran'](vabun.getShp()).to_crs(fun_tiles.proj1),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'vabun_v015',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    if 'albosurv' not in tile_dat:
        print("Loading tiles: albosurv...",flush=True)
        cmap = fun_colors.cmaps['albosurv']
        tile_dat['albosurv'] = {
            'fun': fun_tiles.getShpTiles,
            'dat': cmap['tran'](albosurv.getShp()),
            'cmap': cmap['cmap'],
            'norm': cmap['norm'],
            'label': 'albosurv',
            'cllbl': cmap['cllbl'],
            'clscl': cmap['clscl']
        }
    #
    for prd in papatasi2015.shps:
        if 'papatasi_'+prd not in tile_dat:
            print("Loading tiles: papatasi_%s..." %prd,flush=True)
            cmap = fun_colors.cmaps['papatasi']
            tile_dat['papatasi_'+prd] = {
                'fun': fun_tiles.getShpTiles,
                'dat': cmap['tran'](papatasi2015.shps[prd]).to_crs(fun_tiles.proj1),
                'cmap': cmap['cmap'],
                'norm': cmap['norm'],
                'label': 'papatasi_'+prd,
                'cllbl': cmap['cllbl'],
                'clscl': cmap['clscl']
            }


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

def get_meteo_days(loni, lati, idates, isFeb29):
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

def get_location_ERA5(lon, lat):
    if lon < 0.0:
        lon += 360.0
    #
    loni = fun_clim.getIndex(lon, lwMaskERA5.longitude)
    lati = fun_clim.getIndex(lat, lwMaskERA5.latitude)
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

def get_location(lon, lat, lons, lats):
    if (min(lons) >= 0.0) and (max(lons) >= 180.0) and (lon < 0.0):
        lon += 360.0
    #
    loni = fun_clim.getIndex(lon, lons)
    lati = fun_clim.getIndex(lat, lats)
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

def get_sim_days(loni, lati, idates, isFeb29):
    return {
        "colegg": remove_feb29(annualVectorA.colegg[lati,loni,:],idates,isFeb29),
        "colK": remove_feb29(annualVectorA.colK[lati,loni,:],idates,isFeb29),
        "coln2": remove_feb29(annualVectorA.coln2[lati,loni,:],idates,isFeb29),
        "coln4f": remove_feb29(annualVectorA.coln4f[lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualVectorA.pouts[lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualVectorA.iouts[lati,loni,:],idates,isFeb29)]
    }

def get_sim1980_days(loni, lati, idates, isFeb29):
    return {
        "colegg": remove_feb29(annualVectorA_1980.colegg[lati,loni,:],idates,isFeb29),
        "colK": remove_feb29(annualVectorA_1980.colK[lati,loni,:],idates,isFeb29),
        "coln2": remove_feb29(annualVectorA_1980.coln2[lati,loni,:],idates,isFeb29),
        "coln4f": remove_feb29(annualVectorA_1980.coln4f[lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualVectorA_1980.pouts[lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualVectorA_1980.iouts[lati,loni,:],idates,isFeb29)]
    }

def get_nasa_ssp245_days(loni, lati, idates, isFeb29):
    ssp = 'ssp245'
    return {
        "colegg": remove_feb29(annualNASA.colegg[ssp][lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualNASA.pouts[ssp][lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualNASA.iouts[ssp][lati,loni,:],idates,isFeb29)]
    }

def get_nasa_ssp585_days(loni, lati, idates, isFeb29):
    ssp = 'ssp585'
    return {
        "colegg": remove_feb29(annualNASA.colegg[ssp][lati,loni,:],idates,isFeb29),
        "pouts": [a*100.0 for a in remove_feb29(annualNASA.pouts[ssp][lati,loni,:],idates,isFeb29)],
        "iouts": [a*4000.0 for a in remove_feb29(annualNASA.iouts[ssp][lati,loni,:],idates,isFeb29)]
    }

def get_fcast_days(loni, lati, date0, date1):
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

def get_papatasi_days(loni, lati, idates, isFeb29):
    return {
        "simL": remove_feb29(papatasi2015.simGERI[:,lati,loni],idates,isFeb29),
        "simH": remove_feb29(papatasi2015.simSTENI[:,lati,loni],idates,isFeb29)
    }

def calc_cut(vec,lim,lab):
    if hasattr(vec, '__iter__'):
        return pandas.cut(vec, bins=lim, include_lowest=False, right=False, labels=lab).tolist()
    return lab[numpy.where(numpy.array(lim) > vec)[0][0] - 1]

def get_risk(ret):
    labels = [0, 1, 2]
    larva_limit = [-numpy.inf, 0.07, 3.6, numpy.inf]
    adult_limit = [-numpy.inf, 0.7, 3.6, numpy.inf]
    # larva_limit = [-numpy.inf, 1.0, 50.0, numpy.inf]
    # adult_limit = [-numpy.inf, 10.0, 50.0, numpy.inf]
    pouts_limit = [-numpy.inf, 1.0, 50.0, numpy.inf]
    iouts_limit = [-numpy.inf, 10.0, 500.0, numpy.inf]
    return {
        "larva": calc_cut(ret["coln2"], larva_limit, labels),
        "adult": calc_cut(ret["colegg"], adult_limit, labels),
        "pouts": calc_cut(ret["pouts"], pouts_limit, labels),
        "iouts": calc_cut(ret["iouts"], iouts_limit, labels)
    }

def get_surv(lon,lat):
    return albosurv.getSurv(lon,lat)

def get_decadal(lon, lat, date0, date1=False, ts=False):
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
                               ret['location']['lat'])
    #
    return ret

def get_sandfly(lon, lat, date0, date1=False, ts=False):
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
