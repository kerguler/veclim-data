import json

import veclim_data_server.fun_server as fun_server

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
                    'labels': fun_server.tile_dat[key]['cllbl']
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