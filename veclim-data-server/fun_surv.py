import numpy
import pandas
import json

import geopandas as gpd
import shapely
import shapely.geometry

from scipy.interpolate import interp1d

from environ import DIR_DATA

def approx(a,b,eps=1e-3):
    return (a>=b-eps) & (a<b+eps)

def get_isoweek(a):
    b = int(a.isocalendar()[1])
    return b if b<53 else 1

def dateRange(gr):
    return (numpy.max(gr["date"])-numpy.min(gr["date"])).days

class gridClass_LD:
    def __init__(self) -> None:
        lonlat = json.load(open("%s/clim/lonlat.json" %(DIR_DATA),"r"))
        self.latitude = numpy.array(lonlat['lat'])
        self.longitude = numpy.array(lonlat['lon'])
        self.longitude[self.longitude > 180.0] -= 360.0
        self.lons, self.lats = numpy.meshgrid(self.longitude, self.latitude, indexing='ij')

class presenceAlbopictus:
    def __init__(self) -> None:
        self.shp = [
            {
                'filename': "%s/surveillance/presence/MosquitoAlert_Aedes_albopictus_2023_GAUL" %(DIR_DATA),
                'grabname': 'ADM2_NAME',
                'grabname_alt': 'ADM1_NAME',
                'dataset': "Mosquito Alert (May 2023)",
                'citation': {
                    "publication": "",
                    "database": "GBIF Occurrence Download (10 May 2023)",
                    "url": "http://webserver.mosquitoalert.com/static/tigapublic/spain.html#/en/",
                    "link": "GBIF.org"
                }
            },{
                'filename': "%s/surveillance/presence/AIMSurv_Aedes_albopictus_2019_VectorNetPLG" %(DIR_DATA),
                'grabname': 'GEO_NAME',
                'grabname_alt': 'GEO_NAME',
                'dataset': "AIMSurv - AIM-COST Action (2019)",
                'citation': {
                    "publication": "Miranda Chueca M Á, Barceló Seguí C (2022). AIMSurv <i>Aedes</i> Invasive Mosquito species harmonized surveillance in Europe. AIM-COST Action. Version 2.3. Universitat de les Illes Balears",
                    "database": "GBIF Occurrence Download (09 March 2023)",
                    "url": "https://www.gbif.org/dataset/03269e13-84ae-430f-990e-f11069413e36",
                    "link": "GBIF.org"
                }
            },{
                'filename': "%s/surveillance/presence/VectorNet_Aedes_albopictus_2024May_VectorNetPLG" %(DIR_DATA),
                'grabname': 'GEO_NAME',
                'grabname_alt': 'GEO_NAME',
                'dataset': "ECDC/EFSA Mosquito Maps (May 2024)",
                'citation': {
                    "publication": "ECDC and EFSA. Mosquito maps [internet]. Stockholm: ECDC; 2024",
                    "database": "ECDC/EFSA Mosquito maps (May 2024)",
                    "url": "https://ecdc.europa.eu/en/disease-vectors/surveillance-and-disease-data/mosquito-maps",
                    "link": "Mosquito Maps"
                }
            },{
                'filename': "%s/surveillance/presence/VectorBase_Aedes_albopictus_2024_GAUL" %(DIR_DATA),
                'grabname': 'ADM2_NAME',
                'grabname_alt': 'ADM1_NAME',
                'dataset': "VectorBase (February 2024)",
                'citation': {
                    "publication": "Amos et al., VEuPathDB: the eukaryotic pathogen, vector and host bioinformatics resource center, Nucleic Acids Research, 2021, gkab929",
                    "database": "VectorBase (February 2024)",
                    "url": "https://vectorbase.org/popbio-map/web/?species=Aedes albopictus&view=abnd",
                    "link": "VectorBase"
                }
            },{
                'filename': "%s/surveillance/presence/Kraemer_Aedes_albopictus_2019_GAUL" %(DIR_DATA),
                'grabname': 'ADM2_NAME',
                'grabname_alt': 'ADM1_NAME',
                'dataset': "The Global Compendium (Kraemer et al. 2015)",
                'citation': {
                    "publication": "Kraemer, M., Sinka, M., Duda, K. et al. The global compendium of <i>Aedes aegypti</i> and <i>Ae. albopictus</i> occurrence. Sci Data 2, 150035 (2015)",
                    "database": "The Global Compendium (Kraemer et al. 2015)",
                    "url": "https://doi.org/10.1038/sdata.2015.35",
                    "link": "The paper"
                }
            }
        ]
        self.filename = "%s/surveillance/presence/albopictus_matrix_presence" %(DIR_DATA)
        self.noname = 'Administrative unit not available'
        for shp in self.shp:
            shp['shp'] = gpd.GeoDataFrame.from_file(shp['filename']+".shp", crs='epsg:4326').to_crs("epsg:4326")
            shp['index'] = shp['shp'].sindex
            #
        self.matrix = numpy.load("%s/surveillance/presence/albopictus_matrix_presence.npy" %(DIR_DATA))
        #
    def search(self, lon, lat, res=[0.125,0.125]):
        if lon>180.0:
            lon = lon - 360.0
        #
        shp = shapely.geometry.box(lon-res[0], 
                                   lat-res[1], 
                                   lon+res[0], 
                                   lat+res[1])
        box = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[shp])       
        #
        ret = []
        for shp in self.shp:
            r = {
                'dataset': shp['dataset'],
                'citation': shp['citation'],
                'locations': []
            }
            possible_matches = shp['shp'].iloc[list(shp['index'].query(box.geometry[0]))]
            matches = possible_matches[possible_matches.intersects(box.geometry[0])]
            for j,s in matches.iterrows():
                nm = s[shp['grabname']]
                if nm == self.noname:
                    nm = s[shp['grabname_alt']]
                r['locations'].append(nm)
            if r['locations']:
                ret.append(r)
        #
        return ret
        #
    def getShp(self, res=[0.125,0.125]):
        try:
            gdf = gpd.GeoDataFrame.from_file(self.filename + ".shp", crs='epsg:4326')
        except:
            grid = gridClass_LD()
            geoms = []
            feats = []
            for lon in grid.longitude:
                for lat in grid.latitude:
                    print(lon,lat)
                    tmp = self.search(lon,lat)
                    if len(tmp) == 0:
                        continue
                    geoms.append(shapely.geometry.box(lon-res[0], 
                                                      lat-res[1], 
                                                      lon+res[0], 
                                                      lat+res[1]))
                    feats.append(len(tmp))
            gdf = gpd.GeoDataFrame({
                'mean': feats,
                'geometry': geoms 
            }, crs='epsg:4326')
            gdf = gpd.GeoDataFrame(geometry=[gdf.unary_union], crs=gdf.crs)
            gdf.to_file(self.filename + ".shp")
        #
        gdf['edgecolor'] = 'none'
        gdf['facecolor'] = '#f15a48'
        gdf['label'] = 'Global presence (2024)'
        return gdf

class VectAbundance:
    def __init__(self):
        self.grid = gridClass_LD()
        #
        self.filename = "%s/surveillance/VectAbundance/Vectabundace_v015" %(DIR_DATA)
        vabo = pandas.read_csv(self.filename + ".csv")
        vabo = vabo.loc[
            (vabo["species"] == "albopictus") &
            (~numpy.isnan(vabo["value"])) &
            (vabo["trap_type"] == "ovitrap")
        ]
        vabo["date"] = pandas.to_datetime(vabo["date"])
        #
        self.vb = pandas.DataFrame({
            "date": vabo["date"],
            "year": [a.year for a in vabo["date"]],
            "month": [a.month for a in vabo["date"]],
            "doy": [a.timetuple().tm_yday for a in vabo["date"]],
            "week": [get_isoweek(a) for a in vabo["date"]],
            "lon": vabo["longitude"],
            "lat": vabo["latitude"],
            "gridLon": [self.grid.longitude[numpy.argmin((l-self.grid.longitude)**2)] for l in vabo["longitude"]],
            "gridLat": [self.grid.latitude[numpy.argmin((l-self.grid.latitude)**2)] for l in vabo["latitude"]],
            "samples": vabo["value"]
        })

class AIMsurv:
    def __init__(self):
        self.grid = gridClass_LD()
        #
        self.filename = "%s/surveillance/AIMsurv/occurrence" %(DIR_DATA)
        vabo = pandas.read_csv(self.filename + ".txt", sep="\t", header=0)
        vabo = vabo.loc[
            (vabo["scientificName"]=="Aedes albopictus (Skuse, 1894)") &
            (vabo["hasCoordinate"]) & 
            (~vabo["hasGeospatialIssues"])
        ]
        vabo["date"] = pandas.to_datetime([a.split("/")[1].split("T")[0] for a in vabo["verbatimEventDate"]])
        #
        self.vb = pandas.DataFrame({
            "date": vabo["date"],
            "year": [a.year for a in vabo["date"]],
            "month": [a.month for a in vabo["date"]],
            "doy": [a.timetuple().tm_yday for a in vabo["date"]],
            "week": [get_isoweek(a) for a in vabo["date"]],
            "lon": vabo["decimalLongitude"],
            "lat": vabo["decimalLatitude"],
            "gridLon": [self.grid.longitude[numpy.argmin((l-self.grid.longitude)**2)] for l in vabo["decimalLongitude"]],
            "gridLat": [self.grid.latitude[numpy.argmin((l-self.grid.latitude)**2)] for l in vabo["decimalLatitude"]],
            "samples": vabo["individualCount"]
        })
        #
        self.vb = self.vb.groupby(["gridLon", "gridLat"]).filter(lambda gr: dateRange(gr) > 90)

class VectorBase:
    def __init__(self):
        self.grid = gridClass_LD()
        #
        self.filename = "%s/surveillance/VectorBase/albopictus_VectorBase_2024_02_20" %(DIR_DATA)
        vabo = pandas.read_csv(self.filename + ".csv")
        vabo["date"] = pandas.to_datetime([a.split("/")[0] for a in vabo["Collection date range"]])
        #
        self.vb = pandas.DataFrame({
            "date": vabo["date"],
            "year": [a.year for a in vabo["date"]],
            "month": [a.month for a in vabo["date"]],
            "doy": [a.timetuple().tm_yday for a in vabo["date"]],
            "week": [get_isoweek(a) for a in vabo["date"]],
            "lon": vabo["Longitudes"],
            "lat": vabo["Latitudes"],
            "gridLon": [self.grid.longitude[numpy.argmin((l-self.grid.longitude)**2)] for l in vabo["Longitudes"]],
            "gridLat": [self.grid.latitude[numpy.argmin((l-self.grid.latitude)**2)] for l in vabo["Latitudes"]],
            "samples": vabo["Specimens collected"]
        })
        #
        self.vb = self.vb[self.vb['date'] > pandas.to_datetime('2009-12-31')]
        self.vb = self.vb.groupby(["gridLon", "gridLat"]).filter(lambda gr: dateRange(gr) > 90)

class albosurv:
    def __init__(self):
        self.filename = "%s/surveillance/presence/albosurv_matrix" %(DIR_DATA)
        #
        self.presence = presenceAlbopictus()
        self.vabun = VectAbundance()
        self.aimsurv = AIMsurv()
        self.vbase = VectorBase()
        #
    def _getSurv(self, vb, lon, lat, win=14):
        if lon > 180.0:
            lon -= 360.0
        #
        obs = vb[approx(vb["gridLon"],lon) & approx(vb["gridLat"],lat)]
        if len(obs) == 0:
            return []
        ss = obs.groupby(["week"])["samples"].mean().sort_index().to_frame()
        #
        daily_values = interp1d((ss.index-1)*7, 
                                ss["samples"].values, 
                                kind='linear', 
                                fill_value=None, 
                                bounds_error=False)(numpy.arange(365))
        daily_values[:((ss.index[0]-1)*7)] = None
        daily_values[((ss.index[-1]+1)*7):] = None
        smo = numpy.convolve(daily_values, numpy.ones(win)/win, mode='same')
        #
        return smo
        #
    def _getShp(self, obj, res=[0.125,0.125]):
        try:
            gdf = gpd.GeoDataFrame.from_file(obj.filename + ".shp", crs='epsg:4326')
        except:
            geoms = []
            feats = []
            for gr in obj.vb.groupby(["gridLon","gridLat"]):
                lon, lat = gr[0]
                geoms.append(shapely.geometry.box(lon-res[0], 
                                                  lat-res[1], 
                                                  lon+res[0], 
                                                  lat+res[1]))
                feats.append(gr[1]['samples'].mean(skipna=True))
            gdf = gpd.GeoDataFrame({
                'mean': feats,
                'geometry': geoms 
            }, crs='epsg:4326')       
            gdf.to_file(obj.filename + ".shp")
        #
        return gdf
        #
    def getSurv(self, lon, lat):
        if lon > 180.0:
            lon -= 360.0
        #
        return {
            'vabun': self._getSurv(self.vabun,lon,lat),
            'aimsurv': self._getSurv(self.aimsurv,lon,lat),
            'vbase': self._getSurv(self.vbase,lon,lat)
        }
        #
    def _getShpD(self):
        prs = self.presence.getShp()
        #
        vbn = self._getShp(self.vabun)
        vbn['edgecolor'] = '#1b3958'
        vbn['facecolor'] = 'none'
        vbn['label'] = 'VectAbundance (2010-2022)'
        #
        aim = self._getShp(self.aimsurv)
        aim['edgecolor'] = '#167997'
        aim['facecolor'] = 'none'
        aim['label'] = 'AIMsurv (2020)'
        #
        vbs = self._getShp(self.vbase)
        vbs['edgecolor'] = '#50c0ad'
        vbs['facecolor'] = 'none'
        vbs['label'] = 'VectorBase (2010-2024)'
        #
        return {
            'prs': prs,
            'vbn': vbn,
            'aim': aim,
            'vbs': vbs            
        }

    def getShp(self):
        obj = self._getShpD()
        #
        return gpd.GeoDataFrame(
            pandas.concat([
                obj['prs'], 
                obj['vbn'],
                obj['aim'],
                obj['vbs']
            ], ignore_index=True),
            geometry='geometry', 
            crs=obj['prs'].crs
        )
        #
    def getMatrix(self, res=[0.125,0.125]):
        try:
            mat = numpy.load(self.filename + ".npy")
        except:
            obj = self._getShpD()
            #
            tran = [
                2, # "VectAbundance (2010-2022)"
                3, # "AIMsurv (2020)"
                4, # "VectorBase (2010-2024)"
                1, # "Global presence (2024)"
            ]
            pres = [
                gpd.GeoDataFrame(prs, geometry='geometry', crs=prs.crs)
                for prs in [
                    obj['vbn'],
                    obj['aim'],
                    obj['vbs'],
                    obj['prs'] 
                ]
            ]
            grid = gridClass_LD()
            mat = numpy.ndarray((len(grid.latitude),len(grid.longitude)),dtype=numpy.float64)
            for lati,lat in enumerate(grid.latitude):
                for loni,lon in enumerate(grid.longitude):
                    print("Processing albosurv matrix... %g,%g" %(lat,lon))
                    box = shapely.geometry.box(lon-0.1*res[0], 
                                               lat-0.1*res[1], 
                                               lon+0.1*res[0], 
                                               lat+0.1*res[1])
                    for i,prs in enumerate(pres):
                        cut = prs.intersects(box)
                        if numpy.any(cut):
                            mat[lati,loni] = tran[i]
                            break
            numpy.save(self.filename + ".npy", mat)
        return mat
    