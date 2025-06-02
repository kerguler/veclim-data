import numpy
import pandas
import json

import geopandas as gpd
import shapely
import shapely.geometry

from environ import DIR_DATA

class gridClass_LD:
    def __init__(self) -> None:
        lonlat = json.load(open("%s/clim/lonlat.json" %(DIR_DATA),"r"))
        self.latitude = numpy.array(lonlat['lat'])
        self.longitude = numpy.array(lonlat['lon'])
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
                'filename': "%s/surveillance/presence/VectorBase_Aedes_albopictus_2023_GAUL" %(DIR_DATA),
                'grabname': 'ADM2_NAME',
                'grabname_alt': 'ADM1_NAME',
                'dataset': "VectorBase (March 2023)",
                'citation': {
                    "publication": "Amos et al., VEuPathDB: the eukaryotic pathogen, vector and host bioinformatics resource center, Nucleic Acids Research, 2021, gkab929",
                    "database": "VectorBase (March 2023)",
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

def approx(a,b,eps=1e-3):
    return (a>=b-eps) & (a<b+eps)

def get_isoweek(a):
    b = int(a.isocalendar()[1])
    return b if b<53 else 1

class gridClass_LD:
    def __init__(self) -> None:
        lonlat = json.load(open("%s/clim/lonlat.json" %(DIR_DATA),"r"))
        self.latitude = numpy.array(lonlat['lat'])
        self.longitude = numpy.array(lonlat['lon'])
        self.lons, self.lats = numpy.meshgrid(self.longitude, self.latitude, indexing='ij')

class VectAbundance:
    def __init__(self):
        self.grid = gridClass_LD()
        #
        self.filename = "%s/surveillance/VectAbundance/Vectabundace_v015.csv" %(DIR_DATA)
        vabo = pandas.read_csv(self.filename)
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
        #
    def getSurv(self, lon, lat, dts):
        obs = self.vb[approx(self.vb["gridLon"],lon) & approx(self.vb["gridLat"],lat)]
        if len(obs) == 0:
            return []
        ss = obs.groupby(["week"])["samples"].mean().to_frame()
        obs = pandas.DataFrame({'date': pandas.date_range(dts[0],dts[1])})
        obs["week"] = [get_isoweek(a) for a in obs["date"]]
        obs = obs.merge(ss, on='week', how='left')
        return obs["samples"].values.tolist()
