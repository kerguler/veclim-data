import json

import veclim_data_server.pkg_tiles as pkg_tiles

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

def respondTiles(start_response, kw):
    tile_dat = pkg_tiles.tile_dat
    #
    if not 'date0' in kw:
        return returnResponse(start_response, 'Missing argument: date0')
    date0 = kw['date0']
    #
    if not 'date1' in kw:
        return returnResponse(start_response, 'Missing argument: date1')
    date1 = kw['date1']
    #
    if not 'pr_x' in kw:
        return returnResponse(start_response, 'Missing argument: pr_x')
    pr_x = kw['pr_x']
    #
    if not 'pr_y' in kw:
        return returnResponse(start_response, 'Missing argument: pr_y')
    pr_y = kw['pr_y']
    #
    if not 'pr_z' in kw:
        return returnResponse(start_response, 'Missing argument: pr_z')
    pr_z = kw['pr_z']
    #
    if not 'pr_v' in kw:
        return returnResponse(start_response, 'Missing argument: pr_v')
    pr_v = kw['pr_v']
    #
    v_label = pr_v
    #if ((date0 != None) and (date1 != None)):
    #    v_label += '_dates'
    #    ret = fun_server.load_tiles_dates(v_label,date0,date1)
    #    if ret:
    #        response_body = json.dumps(ret)
    #        return returnResponse(start_response, response_body)
    #
    if ((pr_z == None) or 
        (pr_x == None) or 
        (pr_y == None) or
        (v_label not in tile_dat)):
            print("tile_dat:")
            print(tile_dat)
            ret = {
                key: {
                    'colors': tile_dat[key]['clscl'],
                    'labels': tile_dat[key]['cllbl']
                }
                for key in tile_dat
            }
            response_body = json.dumps(ret)
            return returnResponse(start_response, response_body)                        
    #
    buff = tile_dat[v_label]['fun'](tile_dat[v_label]['dat'], 
                                    pr_z, pr_x, pr_y, 
                                    cmap  = tile_dat[v_label]['cmap'], 
                                    norm  = tile_dat[v_label]['norm'],
                                    label = tile_dat[v_label]['label'])
    response_body = buff
    return returnResponsePNG(start_response, response_body)