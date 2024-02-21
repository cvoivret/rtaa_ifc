from  datetime import datetime

from timezonefinder import TimezoneFinder    
import os
import pandas as pd
import numpy as np

from .sunposition import sunpos

from OCC.Core.gp import gp_Pnt,gp_Dir,gp_Vec,gp_Trsf,gp_Ax1
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_Plane        

from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r
    
from enum import Enum
# define enumeration with auhorised location
class rtaa_dom(Enum):
    RUN=1
    GUY=2
    MAR=3
    GUA=4
    
class rtaa_orientation_sector(Enum):
    NORD=1
    EST=2
    SUD=3
    OUEST=4

# Lat Lon
class rtaa_domcoord(Enum):
    RUN= (-21.114533,55.532062)
    GUY= (3.951795,-53.07823)
    MAR= (14.635541,-61.022814)
    GUA= (16.192207, -61.272382)
    
    @classmethod
    def closer_dom(cls,lon,lat):
        distances=[]
        for lat_lon in rtaa_domcoord:
            distances.append( haversine(lat_lon.value[1],lat_lon.value[0],lon,lat))
        minpos = distances.index(min(distances))
        
        return cls(list(rtaa_domcoord)[minpos])
            

class rtaa_location_orientation:
   
    def __init__(self):
        pass
    
    def set_dom(longitude,latitude):
        self._dom = rtaa_domcoord.RUN.closer_dom(longitude,latitude)
        
         ## datetime to compute shadow
        
        
    # depdending on the dom, load correct dataframes    
    def get_irradiance(self):
        pass
        
    # depending on dom, get critical masked values

class rtaa_sun_data:

    def __init__(self,dom=rtaa_domcoord.RUN,tn_angle=0.0):
        self._dom=dom
        self._tn_angle=tn_angle
        
        absolute_path = os.path.dirname(__file__)
        relative_path = "../data/meteo_rtaa.xlsx"
        full_path = os.path.join(absolute_path, relative_path)
        # depend on the location 
        meteo=pd.read_excel(full_path,sheet_name=self._dom.name)
        
        dt_index=pd.date_range("1/1/2020","12/31/2020",freq='H',inclusive='right')
        self._irradiance=meteo.assign(time=dt_index.values)
        
    def _orientation_sector(self,to_tn):
        orientation=None
        if(abs(to_tn)<=np.pi/4.):
            orientation=rtaa_orientation_sector.NORD
        elif ( to_tn>np.pi/4.) & (to_tn<=3.*np.pi/4.):
            orientation=rtaa_orientation_sector.EST
        elif ( to_tn<-np.pi/4.) & (to_tn>=-3.*np.pi/4.):
            orientation=rtaa_orientation_sector.OUEST
        elif ( abs(to_tn)>3.*np.pi/4.):
            orientation=rtaa_orientation_sector.SUD
        return orientation

    
    def _critical_mask(self,sector=rtaa_orientation_sector.NORD):
        mask=None
        t=self._irradiance.time.dt
        
        if (self._dom is rtaa_domcoord.GUA) | (self._dom is rtaa_domcoord.MAR):
            mask=(t.strftime("%m-%d")<='11-30')*(t.strftime("%m-%d")>='05-01')
        elif  self._dom is rtaa_domcoord.GUY :
            if sector is rtaa_orientation_sector.SUD :
                mask=(t.strftime("%m-%d")<='11-30')*(t.strftime("%m-%d")>='08-01')
            else :
                mask=(t.strftime("%m-%d")<='11-30')*(t.strftime("%m-%d")>='07-01')
        elif self._dom is rtaa_domcoord.RUN :
            if sector is rtaa_orientation_sector.NORD:
                mask=(t.strftime("%m-%d")<='04-30')*(t.strftime("%m-%d")>='02-01')
            else :
                mask=(t.strftime("%m-%d")<='02-28')*(t.strftime("%m-%d")>='01-01')
        else :       
            print( 'Bad configuration of DOM and sector')
                
        return mask
                
    def _critical_irradiance(self,to_tn=0.0):
        
        sector = self._orientation_sector(to_tn)
        
        mask = self._critical_mask(sector)
        
        return self._irradiance[mask]
        
    def sun_data(self,to_tn_angle=0.0):
        
        
        lat=self._dom.value[0]
        lon=self._dom.value[1]
        
        irr= self._critical_irradiance(to_tn_angle).copy()
        dtindex=pd.DatetimeIndex(irr.time)
        #print(dtindex)
        
        tf = TimezoneFinder()
        tz = tf.timezone_at(lng=lon, lat=lat) 
        dr_proj = dtindex.tz_localize(tz)
        dr_proj_utc = dr_proj.tz_convert("UTC")
        
        
        az_vec,zen_vec=sunpos(dr_proj_utc,lat,lon,0)[:2]
        irr["az"]=az_vec
        irr["zen"]=zen_vec       

        origin = gp_Pnt(0.0,0.0,0.0)
        Xaxis = gp_Ax1(origin,gp_Dir(1.0,0.0,0.0))
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
                       
        sun_to_earth_project=[]
        
        for zen,az in zip(zen_vec,az_vec):
            # rotation around X for altitude/elevation setting
            RotX = gp_Trsf()
            RotX.SetRotation(Xaxis,np.deg2rad(90-zen))
            elev_dir = Yaxis.Transformed(RotX)
            
            # rotation around Z axis for azimuth
            # az is generally given clockwise oriented
            RotZ=gp_Trsf()
            RotZ.SetRotation(Zaxis,np.deg2rad(-az) + self._tn_angle)
            sun_axis=elev_dir.Transformed(RotZ)
            
            sun_direction=sun_axis.Direction()
            
            sun_to_earth_project.append(sun_direction.Reversed())
                
        lvec_day=[]
        for v,z in zip(sun_to_earth_project,zen_vec):
            
            if 90-z>0.0:
                #print(v.Coord(),' ',90-z)
                lvec_day.append(v)
        
        daytime = (90.-zen_vec)>=0.0
        irr_crit_day=irr[daytime]
        #print(irr_crit_day)
        
        return lvec_day,irr_crit_day
        
        #return sun_to_earth_project,zen_vec,az_vec


class building_geoloc:
    def __init__(self):
        self._dom = rtaa_domcoord.RUN
        
        self._set_tn_angle(0.0)
        
    def set_from_ifc(self,ifc_file):
        
        repr_context = ifc_file.by_type('IfcGeometricRepresentationContext',False)
        project_repre = repr_context[0]
        true_north = project_repre.TrueNorth
        tn_X,tn_Y= true_north.DirectionRatios
        
        self._set_tn_vec(tn_X,tn_Y)
        
        ifcsite = ifc_file.by_type('IfcSite')[0]
        h,m,s,ms = ifcsite.RefLatitude
        latitude = h+m/60+(s+ms*1e-6)/3600
        h,m,s,ms = ifcsite.RefLongitude
        longitude= h+m/60+(s+ms*1e-6)/3600
        self._set_coord(latitude,longitude)
        
        
    def _set_coord(self,lat,lon):    
        
        self._latitude = lat
        self._longitude= lon
        
        #print("altitude : ",ifcsite.RefElevation)
                
        self._dom = rtaa_domcoord.RUN.closer_dom(self._longitude,self._latitude)
        print("latitude : ",self._latitude)
        print("longitude: ",self._longitude)
        print("dom      : ",self._dom.name)
    
    def _set_tn_angle(self,tn_angle):
        # update angle and vector from angle
        self._tn_angle     = tn_angle
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
                       
        tn_rotation = gp_Trsf()
        tn_rotation.SetRotation(Zaxis,self._tn_angle)
        
        self._tn_vec = gp_Vec(Yaxis.Direction()).Transformed(tn_rotation)
        
        print("True north vector :",self._tn_vec.Coord())
        print("Angle true North : ",self._tn_angle)
        
    def _set_tn_vec(self,tn_X,tn_Y):
        # update angle and vector from vector coord
        self._tn_vec = gp_Vec(tn_X,tn_Y,0.0)
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Xaxis = gp_Ax1(origin,gp_Dir(1.0,0.0,0.0))
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        self._tn_angle =self._tn_vec.Angle(gp_Vec(Yaxis.Direction()))
        
        print("True north vector :",self._tn_vec.Coord())
        print("Angle true North : ",self._tn_angle)
        
    def face_normal_to_tn(self,face):
        """Compute the angle between a normal face and 
        the true north direction

        Parameters
        ----------
        face : Face
            
        Returns
        -------
        to_tn : float
            The angle between the face normal and the true north direction
            
        """        
        srf = BRep_Tool().Surface(face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(face.Orientation()==1):
            face_norm.Reverse()
                        
        Zvec=gp_Vec(0.0,0.0,1.0)
        
        # projection of the face normal on the XY plane
        projected = Zvec.CrossCrossed(gp_Vec(face_norm),Zvec)
        to_tn = projected.AngleWithRef(self._tn_vec,Zvec)
        return to_tn
        
        
