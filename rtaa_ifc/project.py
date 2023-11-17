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
        """Set the true north vector by the value found in ifc file

        Parameters
        ----------
        ifc_file : ifc file object
            the file describing the project
        

        Returns
        -------
        none
            
        """

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
        
        self._tn_angle =self._tn_vec.Angle(gp_Vec(Yaxis.Direction()))
        
        print("Angle true North : ",self._tn_angle)
        
        
    def update_northing_from_angle(self,new_tn):
        """Set the true north vector with a angle value vector

        Parameters
        ----------
        new_tn : float
            Angle between Yaxis and true north
        

        Returns
        -------
        none
            
        """
        self._tn_angle     = new_tn
        
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
        """Set the location of the project, find the time zone

        Parameters
        ----------
        ifc_file : ifc file object
            the file describing the project
        

        Returns
        -------
        none
            
        """
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
            
    
    def face_orientation_angle_tn(self,face):
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
        
        # projection of normal on the XY plane
        
        Zvec=gp_Vec(0.0,0.0,1.0)
        
        projected = Zvec.CrossCrossed(gp_Vec(face_norm),Zvec)
        #print("projected ",projected.Coord())
        to_tn = projected.AngleWithRef(self._tn_vec,Zvec)
        #print(" angle_tn ",to_tn)
        return to_tn
        
    def face_orientation_sector(self,face):
        """Compute the orientation of a face according to the RTAA.
        
        Parameters
        ----------
        face : Face
                    

        Returns
        -------
        orientation : str
            
        """        
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
        """Depending on the orientation of the face, return the table of
        datetime needed to compute the shading factor according to the RTAA

        Parameters
        ----------
        face : Face
            the file describing the project
        Returns
        -------
        mask : np.array(bool)
            
        """        
        
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
       
     
    def data_critical_period_day(self,face):
        """Compute the irradiance and sun position vector for the critical 
        period for the face 

        Parameters
        ----------
        face : Face
            the file describing the project
        

        Returns
        -------
        ( l_vec_day, irr_crit_day): (list(Gp_Vec),pd.DataFrame)
            For day hours : list of sun position vector and irradiance
            
        """    
        
        mask = self.mask_critical_period(face)
        
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
        
        return lvec_day,irr_crit_day
        
     
    
    def sun_vectors(self,dtindex):
        """Compute the sun vector for given datetime. 

        Parameters
        ----------
        dtindex : pandas.DatetimeIndex
            the instants to compute solar position
        

        Returns
        -------
        sun_to_earth_project : list(gp_Vec)
            list of vector
        zen_vec : list(float)
            list of zenith values
        az_vec : list(float)
            list of azimuth values
            
        """
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
