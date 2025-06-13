# http://wsgi.tutorial.codepoint.net/environment-dictionary
# Python's bundled WSGI server
import wsgiserver
from html import escape
from urllib.parse import parse_qs
import json

from environ import VEC_HOST, VEC_PORT

import fun_server

fcast_time = None
def set_fcast_time():
    global fcast_time
    #
    now = fun_server.pandas.Timestamp.today()
    if (fcast_time == None) or ((now - fcast_time).days > 0):
        fun_server.load_forecast_var(reload=True)
        fcast_time = now
        print("Forecast updated",now)

empty_response_PNG = ''
def returnResponsePNG(start_response, response_body):
    status = '200 OK'
    response_headers = [
        ('Access-Control-Allow-Origin', '*'),
        ('Access-Control-Allow-Methods', 'GET'),
        ('Access-Control-Max-Age', '3600'),
        ('Access-Control-Allow-Headers',
         'Content-Type, Content-Length, Access-Control-Allow-Headers, Authorization, X-Requested-With'),
        ('Content-Type', 'image/webp'),
        ('Cache-Control', 'public, max-age=31536000'),
        ('Content-Length', str(len(response_body)))
    ]
    start_response(status, response_headers)
    return [response_body]

empty_response = ''
def returnResponse(start_response, response_body):
    status = '200 OK'
    response_headers = [
        ('Access-Control-Allow-Origin', '*'),
        ('Access-Control-Allow-Methods', 'GET'),
        ('Access-Control-Max-Age', '3600'),
        ('Access-Control-Allow-Headers',
         'Content-Type, Content-Length, Access-Control-Allow-Headers, Authorization, X-Requested-With'),
        ('Content-Type', 'application/json'),
        ('Content-Length', str(len(response_body)))
    ]
    start_response(status, response_headers)
    return [response_body.encode()]

def respondTiles(date0,date1,pr_x,pr_y,pr_z,pr_v,start_response):
    fun_server.load_tiles()
    v_label = pr_v
    if ((date0 != None) and (date1 != None)):
        v_label += '_dates'
        ret = fun_server.load_tiles_dates(v_label,date0,date1)
        if ret:
            response_body = json.dumps(ret)
            return returnResponse(start_response, response_body)
    #
    if ((pr_z == None) or 
        (pr_x == None) or 
        (pr_y == None) or
        (v_label not in fun_server.tile_dat)):
            ret = {
                key: {
                    'colors': fun_server.tile_dat[key]['clscl'],
                    'labels': fun_server.tile_dat[key]['cllbl'],
                    'continuous': fun_server.tile_dat[key]['cont']
                }
                for key in fun_server.tile_dat
            }
            response_body = json.dumps(ret)
            return returnResponse(start_response, response_body)                        
    #
    buff = fun_server.tile_dat[v_label]['fun'](fun_server.tile_dat[v_label]['dat'], 
                                            pr_z, pr_x, pr_y, 
                                            cmap=fun_server.tile_dat[v_label]['cmap'], 
                                            norm=fun_server.tile_dat[v_label]['norm'],
                                            label=fun_server.tile_dat[v_label]['label'])
    response_body = buff
    return returnResponsePNG(start_response, response_body)

def respondAlbopictus(date0,date1,lon,lat,timeseries,meteo_key,sim_key,fcast_key,risk_key,start_response):
    simclm = fun_server.get_decadal(lon,lat,date0,date1,ts=timeseries)
    if not simclm:
        return returnResponse(start_response, empty_response)
    #
    vec = fun_server.presence.search(lon,lat)
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

def respondPapatasi(date0,date1,lon,lat,timeseries,sim_key,start_response):
    simclm = fun_server.get_sandfly(lon,lat,date0,date1,ts=timeseries)
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

def application(environ, start_response):
    # Receive the request from the client (method = GET)
    # ---------------------------------------------------

    parameters = parse_qs(environ.get('QUERY_STRING', ''))
    #
    if 'v' in parameters:
        pr_v = escape(parameters.get('v', [''])[0])
    else:
        pr_v = None
    #
    if 'z' in parameters:
        pr_z = escape(parameters.get('z', [''])[0])
        pr_z = int(pr_z)
    else:
        pr_z = None
    #
    if 'x' in parameters:
        pr_x = escape(parameters.get('x', [''])[0])
        pr_x = int(pr_x)
    else:
        pr_x = None
    #
    if 'y' in parameters:
        pr_y = escape(parameters.get('y', [''])[0])
        pr_y = int(pr_y)
    else:
        pr_y = None
    #
    if 'date' in parameters:
        date0 = escape(parameters.get('date', [''])[0])
        date0 = fun_server.fun_clim.datetime.strptime(date0, "%Y-%m-%d")
        date1 = date0
    elif 'dates' in parameters:
        dates = escape(parameters.get('dates', [''])[0]).split(':')
        date0 = fun_server.fun_clim.datetime.strptime(dates[0], "%Y-%m-%d")
        date1 = fun_server.fun_clim.datetime.strptime(dates[1], "%Y-%m-%d")
    else:
        date0 = None
        date1 = None
        dates = None
    #
    if pr_v != None:
        return respondTiles(date0,date1,pr_x,pr_y,pr_z,pr_v,start_response)
    #
    if 'vec' in parameters:
        vector = escape(parameters.get('vec', [''])[0])
    else:
        vector = 'albopictus'
    if not (vector in fun_server.veclist):
        return returnResponse(start_response, empty_response)
    #
    if 'lon' in parameters:
        lon = escape(parameters.get('lon', [''])[0])
        lon = float(lon)
    else:
        lon = None
    #
    if 'lat' in parameters:
        lat = escape(parameters.get('lat', [''])[0])
        lat = float(lat)
    else:
        lat = None
    #
    if ((lon == None) or 
        (lat == None) or
        (('date' == None) and ('dates' == None))):
        return returnResponse(start_response, empty_response)
    #
    timeseries = False
    meteo_key = 'meteo-mean'
    sim_key = 'sim-mean'
    fcast_key = 'fcast-mean'
    risk_key = 'risk-mean'
    if 'opr' in parameters:
        opr = escape(parameters.get('opr', [''])[0])
        if opr == 'ts':
            timeseries = True
            meteo_key = 'meteo-ts'
            sim_key = 'sim-ts'
            fcast_key = 'fcast-ts'
            risk_key = 'risk-ts'
        else:
            opr = ''
    #
    # Check if forecast should be refreshed
    # ------------------------------------

    set_fcast_time()

    #
    # Process request and respond properly
    # ------------------------------------

    if vector == 'albopictus':
        return respondAlbopictus(date0,date1,lon,lat,timeseries,meteo_key,sim_key,fcast_key,risk_key,start_response)
        #
    elif vector == 'papatasi':
        return respondPapatasi(date0,date1,lon,lat,timeseries,sim_key,start_response)

    return returnResponse(start_response, empty_response)


# Instantiate the server (add certfile and keyfile for SSL)
httpd = wsgiserver.WSGIServer(application,
                              host=VEC_HOST,
                              port=int(VEC_PORT))
#
print("Preparing the datasets...",flush=True)
fun_server.load_global_var()
fun_server.load_tiles()
set_fcast_time()
print("Ready to receive a request...",flush=True)
httpd.start()
