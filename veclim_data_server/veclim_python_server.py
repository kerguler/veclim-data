# http://wsgi.tutorial.codepoint.net/environment-dictionary
# Python's bundled WSGI server
import wsgiserver
from html import escape
from urllib.parse import parse_qs
import json
import pandas

from datetime import datetime

from veclim_data_server.environ import VEC_HOST, VEC_PORT

import veclim_data_server.fun_server as fun_server

import veclim_data_server.response as response

fcast_time = None
def set_fcast_time():
    global fcast_time
    #
    now = pandas.Timestamp.today()
    if (fcast_time == None) or ((now - fcast_time).days > 0):
        fun_server.pkg_sims.reload_forecast_var()
        fcast_time = now
        print("Forecast updated",now)

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
        date0 = datetime.strptime(date0, "%Y-%m-%d")
        date1 = date0
    elif 'dates' in parameters:
        dates = escape(parameters.get('dates', [''])[0]).split(':')
        date0 = datetime.strptime(dates[0], "%Y-%m-%d")
        date1 = datetime.strptime(dates[1], "%Y-%m-%d")
    else:
        date0 = None
        date1 = None
        dates = None
    #
    if 'vec' in parameters:
        vector = escape(parameters.get('vec', [''])[0])
    else:
        vector = 'albopictus'
    if not (vector in fun_server.veclist):
        return response.returnResponse(start_response, response.empty_response)
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
        return response.returnResponse(start_response, response.empty_response)
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

    kw = {
        'date0'         : date0,
        'date1'         : date1,
        'lon'           : lon,
        'lat'           : lat,
        'timeseries'    : timeseries,
        'meteo_key'     : meteo_key,
        'sim_key'       : sim_key,
        'fcast_key'     : fcast_key,
        'risk_key'      : risk_key,
        'start_response': start_response
    }

    if vector in fun_server.pkg_models.modules:
        try:
            return fun_server.pkg_models.modules[vector].respond(**kw)
        except:
            print("ERROR: Problem encountered with request to %s:" %vector, flush=True)
            print(kw, flush=True)

    return response.returnResponse(start_response, response.empty_response)


# Instantiate the server (add certfile and keyfile for SSL)
httpd = wsgiserver.WSGIServer(application,
                              host=VEC_HOST,
                              port=int(VEC_PORT))
#
print("Preparing the datasets...",flush=True)
set_fcast_time()
print("Ready to receive a request...",flush=True)
httpd.start()
