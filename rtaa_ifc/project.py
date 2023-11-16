from  datetime import datetime

from timezonefinder import TimezoneFinder    
import os
import pandas as pd
import numpy as np
import sunposition as sunpos

from OCC.Core.gp import gp_Pnt,gp_Dir,gp_Vec,gp_Pln,gp_Lin,gp_Trsf,gp_Ax3,gp_Ax1
from OCC.Core.gp import gp_Pnt2d,gp_Dir2d,gp_Vec2d,gp_Lin2d

from OCC.Core.BRep import BRep_Tool

from OCC.Core.Geom import Geom_Plane,Geom_Line        
        
class project_location:
    def __init__(self):
        pass
        
    def set_northing_from_ifc(self,ifc_file):
        repr_context = ifc_file.by_type('IfcGeometricRepresentationContext',False)
        project_repre = repr_context[0]
        true_north = project_repre.TrueNorth
        tn_X,tn_Y= true_north.DirectionRatios
        
        # true north vector in project coordinate system
        self._tn_vec = gp_Vec(tn_X,tn_Y,0.0)
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Xaxis = gp_Ax1(origin,gp_Dir(1.0,0.0,0.0))
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        # transformation to apply to convert in project coordinates
        # any vector expressed in world coordinate (sun direction)
        #self._tn_angle_sgn = self._tn_vec.AngleWithRef(gp_Vec(Yaxis.Direction()),gp_Vec(Zaxis.Direction()))
        self._tn_angle =self._tn_vec.Angle(gp_Vec(Yaxis.Direction()))
        
        print("Angle true North : ",self._tn_angle)
        
        #print("Signed angle : ", self._tn_angle_sgn)
        """
        self._tn_rotation = gp_Trsf()
        self._tn_rotation.SetRotation(Zaxis,self._tn_angle)
        print("TN vector : ",(Yaxis.Direction().Transformed(self._tn_rotation).Coord()))
        """
        
    def update_northing_from_angle(self,new_tn):
        #self._tn_angle_sgn = new_tn
        self._tn_angle     = np.abs(new_tn)
        
        print("Angle true North (updated): ",self._tn_angle)
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        #print("
        
        tn_rotation = gp_Trsf()
        tn_rotation.SetRotation(Zaxis,self._tn_angle)
        self._tn_vec = gp_Vec(Yaxis.Direction()).Transformed(tn_rotation)
        print("Updated angle true North : ",self._tn_vec)
        #print("Updated Signed angle : ", self._tn_angle_sgn)
        
        
    
    def set_location_from_ifc(self,ifc_file):
        ## Site location
        ifcsite = ifc_file.by_type('IfcSite')[0]
        h,m,s,ms = ifcsite.RefLatitude
        self._latitude = h+m/60+(s+ms*1e-6)/3600
        h,m,s,ms = ifcsite.RefLongitude
        self._longitude= h+m/60+(s+ms*1e-6)/3600
        print("latitude : ",self._latitude)
        print("longitude: ",self._longitude)
        print("altitude : ",ifcsite.RefElevation)
        
        self._region = 'reunion'
        
        ## datetime to compute shadow
        tf=tf = TimezoneFinder()
        self._tz = tf.timezone_at(lng=self._longitude, lat=self._latitude) 
        print("TimeZone : ",self._tz)

    def infer_rtaa_region_from_ifc(self,ifcfile):
        # set region based on latitude longitude of the project
        pass
        
    def set_rtaa_region(self,region_name='reunion'):
        # test against a list of name
        pass
    
    def get_critical_period_mask(self):
        pass
    
    def load_irradiance(self):
        absolute_path = os.path.dirname(__file__)
        relative_path = "../data/meteo_rtaa.xlsx"
        full_path = os.path.join(absolute_path, relative_path)
        # depend on the location 
        meteo=pd.read_excel(full_path)
        
        #self._dt_index=pd.date_range("1/1/2020","12/31/2020",freq='H',inclusive='right')
        dt_index=pd.date_range("1/1/2020","12/31/2020",freq='H',inclusive='right')

        self._irradiance=meteo.assign(time=dt_index.values)
        
        
        """
        self._critical_period_mask={}
        self._critical_period_mask['reunion']={
        'north':(self._dt_index.strftime("%m-%d")<='04-30')*(self._dt_index.strftime("%m-%d")>='02-01'),
        'ew':(self._dt_index.strftime("%m-%d")<='02-28')*(self._dt_index.strftime("%m-%d")>='01-01'),
        'south':(self._dt_index.strftime("%m-%d")<='02-28')*(self._dt_index.strftime("%m-%d")>='01-01')
        }
        self._critical_period_mask['guyane']={
        'north':(self._dt_index.strftime("%m-%d")<='11-30')*(self._dt_index.strftime("%m-%d")>='07-01'),
        'ew':(self._dt_index.strftime("%m-%d")<='11-30')*(self._dt_index.strftime("%m-%d")>='07-01'),
        'south':(self._dt_index.strftime("%m-%d")<='11-30')*(self._dt_index.strftime("%m-%d")>='08-01')
        }
        """
    """        
    def sunpos(self,dt_index):
        print(self._irradiance)
        vectors, zen_vec, az_vec = self.sun_vectors(dt_index)
        irradiance=self._irradiance[dt_index]
        print(self._dt_index)
        self._irradiance['zenith']=zen_vec
        self._irradiance['azimuth']=az_vec
        
        return vectors, self._irradiance
    """           
    
    def face_orientation_angle_tn(self,face):
        srf = BRep_Tool().Surface(face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(face.Orientation()==1):
            face_norm.Reverse()
        
        # projection of normal on the XY plane
        
        Zvec=gp_Vec(0.0,0.0,1.0)
        
        projected = Zvec.CrossCrossed(gp_Vec(face_norm),Zvec)
        #print("projected ",projected.Coord())
        to_tn = projected.AngleWithRef(self._tn_vec,Zvec)
        #print(" angle_tn ",to_tn)
        return to_tn
        
    def face_orientation_sector(self,face):
        to_tn=self.face_orientation_angle_tn(face)
        
        orientation=None
        if(abs(to_tn)<=np.pi/4.):
            orientation='nord'
        elif ( to_tn>np.pi/4.) & (to_tn<=3.*np.pi/4.):
            orientation='est'
        elif ( to_tn<-np.pi/4.) & (to_tn>=-3.*np.pi/4.):
            orientation='ouest'
        elif ( abs(to_tn)>3.*np.pi/4.):
            orientation='sud'
        return orientation
        
    
    def mask_critical_period(self,face):
        
        # as a function of region return the correct list of sun_position
        # compute face normal and determine the time mask
        to_tn=self.face_orientation_angle_tn(face)
        
        orientation=None
        if(abs(to_tn)<=np.pi/4.):
            orientation='north'
        elif ( (abs(to_tn)>np.pi/4.) & (abs(to_tn)<=3.*np.pi/4.)):
            orientation='ew'
        elif ( abs(to_tn)>3.*np.pi/4.):
            orientation='south'
            
            
        print(' orientation mask ',orientation )
        t=self._irradiance.time.dt
        
        #print(t)
        if self._region == 'reunion':
            if orientation == 'north':
                return (t.strftime("%m-%d")<='04-30')*(t.strftime("%m-%d")>='02-01')
            elif orientation == 'ew':
                return (t.strftime("%m-%d")<='02-28')*(t.strftime("%m-%d")>='01-01')
            elif orientation == 'south' :
                return (t.strftime("%m-%d")<='02-28')*(t.strftime("%m-%d")>='01-01')
        else :
            print("NOT implemented ")
            return
        """
        self._critical_period_mask={}
        self._critical_period_mask['reunion']={
        'north':(self._dt_index.strftime("%m-%d")<='04-30')*(self._dt_index.strftime("%m-%d")>='02-01'),
        'ew':(self._dt_index.strftime("%m-%d")<='02-28')*(self._dt_index.strftime("%m-%d")>='01-01'),
        'south':(self._dt_index.strftime("%m-%d")<='02-28')*(self._dt_index.strftime("%m-%d")>='01-01')
        }
        self._critical_period_mask['guyane']={
        'north':(self._dt_index.strftime("%m-%d")<='11-30')*(self._dt_index.strftime("%m-%d")>='07-01'),
        'ew':(self._dt_index.strftime("%m-%d")<='11-30')*(self._dt_index.strftime("%m-%d")>='07-01'),
        'south':(self._dt_index.strftime("%m-%d")<='11-30')*(self._dt_index.strftime("%m-%d")>='08-01')
        }
                
        return self._critical_period_mask['reunion'][orientation]
        """
        
        """
        # secteur angulaire
        mask=self._critical_period_mask['reunion'][orientation]
        # filtering
        v,_=self.irr_sunpos()
        mask_v = list(np.argwhere(mask).flatten())
        #print(v)
        return list(compress(v,mask_v)),self._irradiance[mask]
        """
     
    def data_critical_period_day(self,face):
    
        
        mask = self.mask_critical_period(face)
        """
        dt_index=pd.DatetimeIndex(self._irradiance.time)
        lvec,zen_vec,az_vec= self.sun_vectors(dt_index)
        self._irradiance['az']=az_vec
        self._irradiance['zen']=zen_vec
        
        for v,z in zip(lvec,zen_vec):
            print(v.Coord(),' ',90-z)
        daytime = (90.-zen_vec)>=0.0
        """
        
        
        
        #dt = self._irradiance[mask]
        irr_crit = self._irradiance[mask].copy()
        
        #print(self._irradiance)
        #print(irr_crit)
        dt_index=pd.DatetimeIndex(irr_crit.time)
        lvec,zen_vec,az_vec= self.sun_vectors(dt_index)
        
        irr_crit["az"]=az_vec
        irr_crit["zen"]=zen_vec
        
        lvec_day=[]
        for v,z in zip(lvec,zen_vec):
            
            if 90-z>0.0:
                #print(v.Coord(),' ',90-z)
                lvec_day.append(v)
        
        daytime = (90.-zen_vec)>=0.0
        irr_crit_day=irr_crit[daytime]
        #print(irr_crit_day)
        """
        # filtering vectors of sun direction
        #mask_v = list(np.argwhere(daytime).flatten())
        print( 'daytime shape ',daytime.shape)
        print(' len mask v ',len(mask_v))
        #lvec_day=list(compress(lvec,mask_v))
        print( 'lvecday  shape ',len(lvec_day))
        print(' zen vec day ',zen_vec[daytime].shape)
        
        for v,z in zip(lvec_day):
            print(v.Coord(),' ',90-z)
        # filtering irradiance
        irr_crit_day=irr_crit[daytime]
        
        print(irr_crit_day)
        """
        #cdcs
        return lvec_day,irr_crit_day
        
        """
        v,irradiance = self.data_critical_period(face)
        print(" critical period date number ",irradiance.shape[0])
        daytime= (90.-irradiance['zenith'])>=0.0
        print(" critical period daylight number ",irradiance[daytime].shape[0])
        #print(daytime[:20])
        mask_v = list(np.argwhere(daytime.values).flatten())
        return list(compress(v,mask_v)),irradiance[daytime]
        """
    
    def sun_vectors(self,dtindex):
        # compute the project sun position from a given time serie (local time, without TZ)
        dr_proj = dtindex.tz_localize(self._tz)
        dr_proj_utc = dr_proj.tz_convert("UTC")
        
        az_vec,zen_vec=sunpos.sunpos(dr_proj_utc,self._latitude,self._longitude,0)[:2]
        
        #elev_vec=90-zen_vec
        
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
        
        
        return sun_to_earth_project,zen_vec,az_vec
