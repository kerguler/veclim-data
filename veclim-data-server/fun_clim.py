import numpy
import xarray
import json
import re
import os
from datetime import datetime, timedelta, date

import geopandas as gpd
import shapely
import shapely.geometry

from environ import DIR_DATA

def stdwarn(msg):
    print("ERROR",msg,flush=True)

def is_leap_year(y):
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False

def getIndex(vec,v):
    return numpy.argmin(numpy.abs(vec-v))

def daylength(lat,day):
    """
    Translated from the daylength function of the geosphere package of R
    lat: latitude in degree decimal (float)
    day: datetime.date or day of the year (integer)
    """
    if isinstance(day,date):
        day = day.timetuple().tm_yday
    pi180 = numpy.pi / 180.0
    P = numpy.arcsin(0.39795 * numpy.cos(0.2163108 + 2 * numpy.arctan(0.9671396 * numpy.tan(0.0086 * (day - 186)))))
    a = (numpy.sin(0.8333 * pi180) + numpy.sin(lat * pi180) * numpy.sin(P))/(numpy.cos(lat * pi180) * numpy.cos(P))
    a = numpy.min([numpy.max([a, -1]), 1])
    DL = 24 - (24.0/numpy.pi) * numpy.arccos(a)
    return(DL)

class lwMaskERA5:
    def __init__(self) -> None:
        self.lwmask = numpy.load("%s/clim/lwmask_0.1lw.npy" %(DIR_DATA))
        lonlat = json.load(open("%s/clim/lonlat.json" %(DIR_DATA),"r"))
        self.latitude = numpy.array(lonlat['lat'])
        self.longitude = numpy.array(lonlat['lon'])
        #
    def island(self,loni,lati):
        return ~self.lwmask[lati,loni]

class annualERA5:
    def __init__(self) -> None:
        self.decade = "2010-2019"
        dr = "%s/clim/ERA5/ERA5_single_levels_decadal/2010_to_2019" %(DIR_DATA)
        atemp = xarray.open_dataset("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_temperature.nc" %(dr))
        tf = 365
        latitude = atemp['lat'].values
        self.photo = numpy.array([[daylength(lat, d) for d in numpy.arange(tf)] for lat in latitude])
        self.atemp = atemp['2m_temperature'].values
        self.atmin = xarray.open_dataset("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_temperature_min.nc" %(dr))['2m_temperature_min'].values
        self.atmax = xarray.open_dataset("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_temperature_max.nc" %(dr))['2m_temperature_max'].values
        self.rehum = xarray.open_dataset("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_2m_relative_humidity.nc" %(dr))['2m_relative_humidity'].values
        self.precp = xarray.open_dataset("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_total_precipitation.nc" %(dr))['total_precipitation'].values
        self.soilw = xarray.open_dataset("%s/ERA5_ERA5_single_levels_decadal_2010_to_2019_volumetric_soil_water_layer_1.nc" %(dr))['volumetric_soil_water_layer_1'].values

class annualVectorA:
    def __init__(self) -> None:
        self.decade = "2010-2019"
        self.dlabel = "2010_to_2019"
        dr = "%s/sims/vector08c_Q4.a100+1/ERA5_single_levels_decadal/%s" %(DIR_DATA,self.dlabel)
        self.colegg = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_colegg.nc" %(dr,self.dlabel))['colegg'].values
        self.colK = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_colK.nc" %(dr,self.dlabel))['colK'].values
        self.coln2 = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_coln2.nc" %(dr,self.dlabel))['coln2'].values
        self.coln4f = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_coln4f.nc" %(dr,self.dlabel))['coln4f'].values
        #
        dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI/ERA5_single_levels_decadal/p4000r100w60/%s" %(DIR_DATA,self.dlabel)
        self.iouts = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_chikv_QI_ERA5_single_levels_decadal_p4000r100w60_%s_iouts.nc" %(dr,self.dlabel))['iouts'].values
        self.pouts = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_chikv_QI_ERA5_single_levels_decadal_p4000r100w60_%s_pouts.nc" %(dr,self.dlabel))['pouts'].values
        #
        tmp = numpy.genfromtxt("%s/surveillance/Italy2008/coord_albopictus_Italy2008.csv" %(DIR_DATA), delimiter=',',names=True)
        self.acc = {}
        for row in tmp:
            loni = int(row['loni'])
            lati = int(row['lati'])
            if not lati in self.acc:
                self.acc[lati] = {}
            self.acc[lati][loni] = 1
        #
    def accuracy(self,lati,loni):
        if not lati in self.acc or not loni in self.acc[lati]:
            return 0
        return self.acc[lati][loni]

class annualVectorA_1980:
    def __init__(self) -> None:
        self.decade = "1980-1989"
        self.dlabel = "1980_to_1989"
        #
        dr = "%s/sims/vector08c_Q4.a100+1/ERA5_single_levels_decadal/%s" %(DIR_DATA,self.dlabel)
        self.colegg = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_colegg.nc" %(dr,self.dlabel))['colegg'].values
        self.colK = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_colK.nc" %(dr,self.dlabel))['colK'].values
        self.coln2 = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_coln2.nc" %(dr,self.dlabel))['coln2'].values
        self.coln4f = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_ERA5_single_levels_decadal_%s_coln4f.nc" %(dr,self.dlabel))['coln4f'].values
        #
        dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI/ERA5_single_levels_decadal/p4000r100w60/%s" %(DIR_DATA,self.dlabel)
        self.iouts = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_chikv_QI_ERA5_single_levels_decadal_p4000r100w60_%s_iouts.nc" %(dr,self.dlabel))['iouts'].values
        self.pouts = xarray.open_dataset("%s/sims_vector08c_Q4.a100+1_chikv_QI_ERA5_single_levels_decadal_p4000r100w60_%s_pouts.nc" %(dr,self.dlabel))['pouts'].values

class annualNASA:
    def __init__(self) -> None:
        self.models = ['ACCESS-CM2','ACCESS-ESM1-5','EC-Earth3']
        self.ssps = ['ssp245','ssp585']
        self.decade = "2090-2099"
        self.dlabel = "2090_to_2099"
        dr = "%s/sims/vector08c_Q4.a100+1" %(DIR_DATA)
        self.colegg = {
            ssp: numpy.nanmean([
                xarray.open_dataset("%s/NASA_decadal_%s_%s/%s/sims_vector08c_Q4.a100+1_NASA_decadal_%s_%s_%s_colegg.nc" %(dr,model,ssp,self.dlabel,model,ssp,self.dlabel))['colegg'].values
                for model in self.models
            ], axis=0)
            for ssp in self.ssps
        }
        dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI" %(DIR_DATA)
        self.iouts = {
            ssp: numpy.nanmean([
                xarray.open_dataset("%s/NASA_decadal_%s_%s/p4000r100w60/%s/sims_vector08c_Q4.a100+1_chikv_QI_NASA_decadal_%s_%s_p4000r100w60_%s_iouts.nc" %(dr,model,ssp,self.dlabel,model,ssp,self.dlabel))['iouts'].values
                for model in self.models
            ], axis=0)
            for ssp in self.ssps
        }
        self.pouts = {
            ssp: numpy.nanmean([
                xarray.open_dataset("%s/NASA_decadal_%s_%s/p4000r100w60/%s/sims_vector08c_Q4.a100+1_chikv_QI_NASA_decadal_%s_%s_p4000r100w60_%s_pouts.nc" %(dr,model,ssp,self.dlabel,model,ssp,self.dlabel))['pouts'].values
                for model in self.models
            ], axis=0)
            for ssp in self.ssps
        }
        #
        self.popsize, self.reps, self.win = [4000, 100, 60]

class forecastECMWF:
    def __init__(self) -> None:
        dr = "%s/sims/vector08c_Q4.a100+1/ECMWF" %(DIR_DATA)
        #
        self.fld = sorted(os.listdir(dr))[-1]
        self.overlap = [re.sub(r"\_","-",s) for s in re.findall(r"\d\d\d\d\_\d+\_\d+",self.fld)]
        #
        if (len(self.overlap)!=2):
            stdwarn("Vector forecast directory is misworded! %s" %self.fld)
            return
        #
        self.colegg = xarray.open_dataset("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_colegg.nc" %(dr,self.fld,self.fld))['colegg'].values
        self.colK = xarray.open_dataset("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_colK.nc" %(dr,self.fld,self.fld))['colK'].values
        self.coln2 = xarray.open_dataset("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_coln2.nc" %(dr,self.fld,self.fld))['coln2'].values
        coln4f = xarray.open_dataset("%s/%s/sims_vector08c_Q4.a100+1_ECMWF_%s_coln4f.nc" %(dr,self.fld,self.fld))
        self.coln4f = coln4f['coln4f'].values
        self.dates = coln4f['time'].values
        #
        dr = "%s/sims/vector08c_Q4.a100+1_chikv_QI/ECMWF" %(DIR_DATA)
        #
        fld = sorted(os.listdir(dr))[-1]
        overlap = [re.sub(r"\_","-",s) for s in re.findall(r"\d\d\d\d\_\d+\_\d+",fld)]
        print("Reading forecast...", overlap)
        #
        if (len(overlap)!=2):
            stdwarn("Disease forecast directory is misworded! %s" %fld)
            return
        if (overlap[0] != self.overlap[0]) or (overlap[1] != self.overlap[1]):
            stdwarn("Dates of vector and disease forcasts mismatch! %s != %s or %s != %s" %(overlap[0],self.overlap[0],overlap[1],self.overlap[1]))
            return
        #
        tmp = re.findall(r"p(\d+?)r(\d+?)w(\d+?)",fld)
        if len(tmp)!=1:
            stdwarn("Disease forecast configuration is missing or misworded! %s" %fld)
            return
        #
        self.popsize, self.reps, self.win = tmp[0]
        #
        self.iouts = xarray.open_dataset("%s/%s/sims_vector08c_Q4.a100+1_chikv_QI_ECMWF_%s_iouts.nc" %(dr,fld,fld))['iouts'].values
        self.pouts = xarray.open_dataset("%s/%s/sims_vector08c_Q4.a100+1_chikv_QI_ECMWF_%s_pouts.nc" %(dr,fld,fld))['pouts'].values
        #
    def getOverlap(self):
        return self.overlap
        #
    def getTransConf(self):
        return {
            'popsize': self.popsize, 
            'reps': self.reps, 
            'win': self.win
        }

class papatasi2015:
    def __init__(self) -> None:
        dr = "%s/sims/sand_papatasi/2015" %(DIR_DATA)
        prds = ['aprdec','aprjun','julsep','octdec']
        lonlat = json.load(open("%s/lonlat_sand.json" %(dr),"r"))
        self.latitude = numpy.array(lonlat['lat'])
        self.longitude = numpy.array(lonlat['lon'])
        self.shps = {
            prd: gpd.read_file("%s/intersect_clc_sim_%s.shp" %(dr,prd))
            for prd in prds
        }
        self.simGERI = numpy.load("%s/mech_model_STENI_papatasi_combinedA_colnvAf_posterior_mean_GERI.npy" %(dr))
        self.simSTENI = numpy.load("%s/mech_model_STENI_papatasi_combinedA_colnvAf_posterior_mean_STENI.npy" %(dr))
        #
        self.mask = numpy.isnan(numpy.nanmean(self.simSTENI[90:,:,:],axis=0))
        #
    def island(self,loni,lati):
        return not self.mask[lati,loni]

class popDens:
    def __init__(self) -> None:
        self.pop = numpy.load("%s/clim/SEDAC/gpw_v4_population_density_adjusted_fromHiRes_0.1lwmask_2010_to_2020.npy" %(DIR_DATA))

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
