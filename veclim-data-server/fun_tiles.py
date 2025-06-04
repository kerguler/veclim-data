import matplotlib
matplotlib.use('agg')

import numpy

import matplotlib.pyplot as plt
import cartopy

from imageio.v3 import imwrite, imread
from io import BytesIO

import os

from environ import DIR_DATA

pxbounds = []
imsize = 256
fig_dpi = 100
proj0 = cartopy.crs.PlateCarree(central_longitude=180.0-0.125)
proj0s = cartopy.crs.epsg(3035)
proj1 = cartopy.crs.Mercator.GOOGLE
# proj1 = cartopy.crs.epsg(3857)
extbound = [-180, 180, -90+0.125, 90+0.125]

coords = {}

def initSubmap():
    fig = plt.figure(figsize=(imsize/fig_dpi,imsize/fig_dpi),tight_layout={'pad':0,'w_pad':0,'h_pad':0},frameon=False)
    fig.patch.set_alpha(0)
    ax = fig.add_subplot(111,projection=proj1)
    ax.patch.set_alpha(0)
    ax.set_axis_off()
    #
    # ax.coastlines()
    #
    tmp = ax.get_extent()
    pxb = [tmp[0], tmp[2], tmp[1], tmp[3]]
    #
    return {
        'fig': fig, 
        'ax': ax, 
        'pxb': pxb
    }

def plotCanvas(fig):
    fig.canvas.draw()
    w,h = fig.canvas.get_width_height()
    buff = fig.canvas.print_to_buffer()[0]
    image = numpy.frombuffer(buff, dtype=numpy.uint8)
    image.shape = (w,h,4)
    return {
        'array': image,
        'buffer': buff
    }

def saveTile(label, x, y, z, mat):
    subpath = DIR_DATA
    for pth in ['tiles', label, str(z), str(x)]:
        subpath += "/"+pth
        if not os.path.exists(subpath):
            try:
                os.makedirs(subpath)
            except:
                pass
    #
    file = "%s/%d.webp" %(subpath, y)
    imwrite(file, mat, format='webp')

def loadTile(label, x, y, z):
    subpath = DIR_DATA
    for pth in ['tiles', label, str(z), str(x)]:
        subpath += "/"+pth
    #
    file = "%s/%d.webp" %(subpath, y)
    if os.path.exists(file):
        return imread(file, extension='.webp')
    else:
        return []

def buffArray(mat):
    buff = BytesIO()
    imwrite(buff, mat, format='webp', extension='.webp')
    buff.seek(0)
    return buff.read()

def getTiles(dat, pr_z, pr_x, pr_y, cmap=None, norm=None, label=''):
    if label:
        carray = loadTile(label,pr_x,pr_y,pr_z)
        if len(carray) > 0:
            return buffArray(carray)
    #
    submap = initSubmap()
    #
    npieces = 2**pr_z
    pxbounds = submap['pxb']
    xw = (pxbounds[2] - pxbounds[0])/npieces
    yh = (pxbounds[3] - pxbounds[1])/npieces
    pxbounds = [
        pxbounds[0] + pr_x*xw,
        pxbounds[3] - (pr_y+1)*yh,
        pxbounds[0] + (pr_x+1)*xw,
        pxbounds[3] - pr_y*yh,
    ]
    #
    submap['ax'].set_xlim([pxbounds[0], pxbounds[2]])
    submap['ax'].set_ylim([pxbounds[1], pxbounds[3]])
    #
    # Note: This has to come after setting the axis limits
    # Note: vmin and vmax should be set (not left as None)
    try:
        submap['ax'].imshow(dat,
                            origin="upper",
                            interpolation='none',
                            transform=proj0,
                            cmap=cmap,
                            norm=norm,
                            extent=extbound)
    except Exception as e:
        pass
    #
    canvas = plotCanvas(submap['fig'])
    buffr = buffArray(canvas['array'])
    plt.close(submap['fig'])
    #
    if label:
        saveTile(label,pr_x,pr_y,pr_z,canvas['array'])
    #
    return buffr

# TO DO
# -----
# This needs to be combined with getTiles 
# -----
def getShpTiles(shp, pr_z, pr_x, pr_y, cmap=None, norm=None, label=''):
    if label:
        carray = loadTile(label,pr_x,pr_y,pr_z)
        if len(carray) > 0:
            return buffArray(carray)
    #
    submap = initSubmap()
    #
    npieces = 2**pr_z
    pxbounds = submap['pxb']
    xw = (pxbounds[2] - pxbounds[0])/npieces
    yh = (pxbounds[3] - pxbounds[1])/npieces
    pxbounds = [
        pxbounds[0] + pr_x*xw,
        pxbounds[3] - (pr_y+1)*yh,
        pxbounds[0] + (pr_x+1)*xw,
        pxbounds[3] - pr_y*yh,
    ]
    #
    submap['ax'].set_xlim([pxbounds[0], pxbounds[2]])
    submap['ax'].set_ylim([pxbounds[1], pxbounds[3]])
    #
    # Note: This has to come after setting the axis limits
    # Note: vmin and vmax should be set (not left as None)
    if cmap == None and norm == None:
        shp.plot(edgecolor='black', 
                 facecolor='none',
                 linewidth=2,
                 ax=submap['ax'])
    else:
        shp.plot(column='mean',
                 cmap=cmap,
                 norm=norm,
                 ax=submap['ax'])
    #
    canvas = plotCanvas(submap['fig'])
    buffr = buffArray(canvas['array'])
    plt.close(submap['fig'])
    #
    if label:
        saveTile(label,pr_x,pr_y,pr_z,canvas['array'])
    #
    return buffr
