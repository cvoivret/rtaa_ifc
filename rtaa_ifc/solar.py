
from collections import defaultdict
from collections import Counter

from itertools import compress,tee,combinations,product
import numpy as np
import scipy.integrate as spi
import pandas as pd

import ifcopenshell

from .project import project_location
from .geom import   (
                    ifcelement_as_solid,
                    fuse_listOfShape,
                    shapes_as_solids,
                    get_external_shell,
                    exterior_wall_normal_and_plane,
                    biggestface_along_vector
                    )

from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color,Quantity_TOC_RGB

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism,BRepPrimAPI_MakeCylinder
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties,brepgprop_VolumeProperties
from OCC.Core.BRepBuilderAPI import (
                                    BRepBuilderAPI_Copy,
                                    BRepBuilderAPI_Transform,
                                    BRepBuilderAPI_MakeFace
                                    )
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepTools import breptools_UVBounds

from OCC.Core.GProp import GProp_GProps
from OCC.Core.GeomLProp import GeomLProp_SLProps

from OCC.Core.gp import (
                        gp_Pnt,
                        gp_Dir,
                        gp_Vec,
                        gp_Pln,
                        gp_Lin,
                        gp_Trsf,
                        gp_Ax3,
                        gp_Ax1
                        )


from OCC.Core.Geom import Geom_Plane

from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.TopTools import TopTools_ListOfShape,TopTools_IndexedMapOfShape
from OCC.Core.TopExp import topexp_MapShapes
from OCC.Core.TopAbs import TopAbs_FACE

from OCC.Extend.TopologyUtils import TopologyExplorer
 
from OCC.Core.BOPAlgo import BOPAlgo_BOP,BOPAlgo_Operation
from OCC.Core.BOPAlgo import BOPAlgo_CellsBuilder



def compute_direct_mask_on_face(sun_dir,mybuilding,myface,theface_norm):
    """Compute the shadow geometry generated by my building on myface 
    trough solid clipping method for one solar vector sun_dir.
    

    Parameters
    ----------
    sun_dir : gp_Dir
    mybuilding : Solid
    myface : Face
    theface_norm : Face (could be removed ?)

    Returns
    -------
    Shadow geometry : Face 
        
    """  
    
    building = BRepBuilderAPI_Copy(mybuilding).Shape()
    theface = BRepBuilderAPI_Copy(myface).Shape()
    # face not exposed to the sun
    if theface_norm.Dot(sun_dir)>-1e-5:
        return theface
    
    gpp=GProp_GProps()
    brepgprop_SurfaceProperties(theface,gpp)
    gf_area=gpp.Mass()
    
    ext_vec=gp_Vec(sun_dir)
    ext_vec.Multiply(10.)
    
    # extrusion 
    extrusion1=BRepPrimAPI_MakePrism(theface,-ext_vec,False,True).Shape()
   
    # intersection
    los2 = TopTools_ListOfShape()
    los_avoid = TopTools_ListOfShape()
    los2.Append(building)
    los2.Append(extrusion1)
    
    
    precision = 1e-9
    ntry=1
    haswarning=True
    
    cb=BOPAlgo_CellsBuilder()
    cb.SetNonDestructive(True)
    cb.SetArguments(los2)
    cb.SetFuzzyValue(precision)
    cb.Perform()
    cb.RemoveAllFromResult()
    cb.AddToResult(los2,los_avoid,2)
    
    cb.RemoveInternalBoundaries()
    haswarning =   cb.HasWarnings()
    # redo intersection with more tolerance
    # stop after 50 try
    while( haswarning ):
        precision*=1.2
        ntry+=1
                            
        cb=BOPAlgo_CellsBuilder()
        cb.SetNonDestructive(True)
        cb.SetArguments(los2)
        cb.SetFuzzyValue(precision)
        cb.Perform()
        cb.RemoveAllFromResult()
        cb.AddToResult(los2,los_avoid,2)
        cb.RemoveInternalBoundaries()
        haswarning =   cb.HasWarnings()
        
        if ( ntry>50):
            break
            
    intersection=cb.Shape()
    
    #Selection of face oriented to the sun and
    # with a significant flux 
    intersection_faces=list(TopologyExplorer(intersection).faces())
    if(len(intersection_faces)==0):
        return TopoDS_Face()           
    # flux trough the face of interest
    # we will compare this flux and the fluw trough 
    # each face of the intersetion to filter 
    brepgprop_SurfaceProperties(theface,gpp)
    theface_area=gpp.Mass()
    theface_flux=theface_norm.Dot(sun_dir)*theface_area
    
    larea=[]
    lfaces=[]
          
    for ff in intersection_faces:
        
        adapt=BRepAdaptor_Surface(ff)
        # normal of cylinder faces is not unique
        # dirty trick to extract fraction of cylinder face that 
        # can cas shadows
        if adapt.GetType()==1:
            cyl=adapt.Cylinder()
            umin,umax,vmin,vmax=breptools_UVBounds(ff)
            if vmin<0.0:
                cyl.VReverse()
            
            ax3=cyl.Position()
            vec=gp_Dir(*sun_dir.Coord())
            
            vec.Cross(ax3.Direction())
            newax3=gp_Ax3(ax3.Location(),ax3.Direction(),vec)
            shape=BRepPrimAPI_MakeCylinder(newax3.Ax2(),cyl.Radius()*2,2,3.14).Shape()
            
            com=BRepAlgoAPI_Common(shape,ff)
            com.Build()
            shape=com.Shape()
            #lcyl.append(shape)
            maps=TopTools_IndexedMapOfShape()
            topexp_MapShapes(shape,TopAbs_FACE,maps)
            lfacetokeep=[maps.FindKey(i) for i in range(1,maps.Size()+1)]
            if( len(lfacetokeep)==1):
                ff=lfacetokeep[0]
            else:
                continue
        
        # normal computation
        srf3 = BRep_Tool().Surface(ff)
        umin,umax,vmin,vmax=breptools_UVBounds(ff)
        props=GeomLProp_SLProps(srf3,0.5*(umax-umin),0.5*(vmax-vmin),1,0.001)
        fn=props.Normal()
        if(ff.Orientation()==1):
            fn.Reverse()
        
        # avoid face nearly parallel with extrusion generatrix
        # ie face with normal perpendicular with extrusion direction
        if(fn.Dot(sun_dir)<-1e-5):
            brepgprop_SurfaceProperties(ff,gpp)
            larea.append(gpp.Mass())
            faceflux = larea[-1]*fn.Dot(sun_dir)
           
            # filtering out faces correctingly oriented but with small flux
            if (faceflux/theface_flux>1e-4):
                lfaces.append(ff)
    
    # extrusion of selected faces
    lsolid=[ BRepPrimAPI_MakePrism(s,ext_vec,False,True).Shape() for s in lfaces]
    
    # void face with zero area
    if(len(lsolid)==0):
        return TopoDS_Face() 
        
    # intersectionS with the face of interest
    lface2=[]
    for s,f in zip(lsolid,lfaces):
        common=BRepAlgoAPI_Common(s,theface)
        sh=common.Shape()
        lface2.append(sh)
    
    
    
    if len(lface2)==1:
        shadowface=lface2[0]
    else:    
        los2 = TopTools_ListOfShape()
        for s in lface2:
            los2.Append(s)
        # Union of faces to obtain one face
        # Cell builder handles well the holes
        cb=BOPAlgo_CellsBuilder()
        cb.SetArguments(los2)
        cb.Perform()
        cb.AddAllToResult(2,False)
        cb.RemoveInternalBoundaries()
        shadowface=cb.Shape()
    
    if shadowface==None:
        shadowface=TopoDS_Face()
    
    
    return shadowface


def compute_diffuse_mask_on_face(mybuilding,theface,theface_norm):
    """

    Parameters
    ----------
    mybuilding : Shape
    theface : Face
    theface_norm : gp_Vec

    
    Returns
    -------
    areas,mask_vector,projected_vec,mask_sky,lvec
        
    """        

    
    # set a discretization of the half sphere exposed to the sun
    n=10
    arct=np.arccos(np.linspace(0,1,n)[::-1])
    t=arct
    
    temp=np.arccos(np.linspace(1,-1,n))
    arcp=np.concatenate((temp,temp[1:]+np.pi))
    p=arcp
    
    # integration points
    mid_t=.5*(t[1:]+t[:-1])
    mid_p=.5*(p[1:]+p[:-1])
    # coordinate to set the direction of mask calculation
    t1,p1=np.meshgrid(mid_t,mid_p,indexing='ij')
    # areas of each solid angle quadrant
    areas = np.zeros_like(t1)
    
    # computation of local area associated with each direction
    # sum = 2*pi
    
    # funtion to integrate to compute the area of each quadrant
    ds=lambda x,y: np.sin(x)
    
    for i,j in product(range(t1.shape[0]),range(t1.shape[1])):    
        res,err = spi.nquad(ds, [[t[i], t[i+1]],[p[j], p[j+1]]])
        areas[i,j]=res
    areas=areas.flatten()   
    
    # To Be Improved :
    # presently the central direction of the dome is set along the 
    # face normal.
    # It is different of many shading factor implementation and
    # I (CV) am not sure if it is pertinent. 
    # It is quite logical but must be check if it is coherent 
    # with the terme (1+beta) for the diffuse part.
    # It work for the RTAAdom becaus of the ratio but it is 
    # suspicious for a plain shading factor.
    #The following rotations could be then avoided
    
    # first incline (between 0+delta (normal) and PI/2) with respect of one of the axis of the surface 
    # we dont care which axis because of the next rotation
    srf = BRep_Tool().Surface(theface)
    adapt=BRepAdaptor_Surface(theface)
    plane = adapt.Plane()
    Xaxis_plane_ax=plane.XAxis()
    normal_plane_ax = plane.Axis()
    normal_plane_vec = gp_Vec(normal_plane_ax.Direction()) # unit vector
    
    origin = gp_Pnt(0.0,0.0,0.0)
    Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
    
    t_transfo = gp_Trsf()
    p_transfo = gp_Trsf()
    
    lvec=[]
    lvec_proj=[]
    lmask_sky=[]
    for rot_t,rot_p in zip(t1.flat,p1.flat):
        t_transfo.SetRotation(Xaxis_plane_ax,rot_t)
        p_transfo.SetRotation(normal_plane_ax,rot_p)
        transformed = theface_norm.Transformed(t_transfo)
        transformed = transformed.Transformed(p_transfo)
        final = gp_Vec(transformed) # oriented from building to sun
        lvec.append(gp_Dir(final.Reversed())) # set direction : from sun to building
        lmask_sky.append( final.Dot(gp_Vec(Zaxis.Direction()))>=0.)
        # projection of the direction on the base plane of the dome
        # by double vectorial product
        projected = normal_plane_vec.CrossCrossed(final,normal_plane_vec)
        lvec_proj.append(projected.Coord())
    
    # remove the axis with null coordinate (in surface plan)
    # to obtain the projection in the base plane coordinate 
    projected_vec= np.array(lvec_proj)
    idx=np.argwhere(np.all(np.isclose(projected_vec[...,:],np.zeros(projected_vec.shape)),axis=0))
    projected_vec=np.delete(projected_vec,idx,axis=1)
    
    
    # compute the direct shading factor for each direction
    dmof=direct_mask_on_face(theface,lvec)
    dmof.compute_mask(mybuilding)
    dmof.shadow_areas()
    dmof.face_area()
    dmof.mask_ratio()
    
    mask_vector = np.array(dmof._mask_ratio)
    mask_sky = np.array(lmask_sky)
    
    return areas,mask_vector,projected_vec,mask_sky,lvec 
    
    
class diffuse_mask_on_face:
    def __init__(self,face):
        self._face=face
    
    def compute_mask(self,exposed_building):
        """Compute the value of shading factor for a given
        discretization of the sky halfdome 

        Parameters
        ----------
        exposed_building : Shape
        min_area : float
            ratio of flux to filter out face in shading 
            factor computation
        Returns
        -------
        The cutted building : Solid
            
        """        
        srf = BRep_Tool().Surface(self._face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
            
              
        areas,mask_vector,vec,sky,lvec=compute_diffuse_mask_on_face(exposed_building,self._face,face_norm)
        
        self._mask_value=mask_vector
        self._masks_areas=areas
        self._directions=vec 
        self._mask_sky=sky
        self._lvec = lvec
        self._ldot = -1.*np.array([ face_norm.Dot(v) for v in lvec])
        
        

    def compute_weighted_mask(self):
        """

        Parameters
        ----------
        None
        Returns
        -------
        None
            
        """        
        
        m=self._mask_value
        a=self._masks_areas
        ms=self._mask_sky
        dot=self._ldot
        
        self._wm=(m*dot*a).sum()/(a*dot).sum()
        self._wm_sky=(m[ms]*a[ms]*dot[ms]).sum()/((a[ms]*dot[ms]).sum())
        self._wm_soil=(m[~ms]*a[~ms]*dot[~ms]).sum()/((a[~ms]*dot[~ms]).sum())

       
        print("\n wm ", self._wm)
        print(" sky ",self._wm_sky)
        print(" soil ", self._wm_soil)
        
        
    def display(self,exposed_building):
        """Display shapes

        Parameters
        ----------
        exposed_building : Shape
        
        Returns
        -------
        None
            
        """        
        
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)

        srf = BRep_Tool().Surface(self._face)
        umin,umax,vmin,vmax=breptools_UVBounds(self._face)
        center = srf.Value((umax+umin)*.5,(vmax+vmin)*.5)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
        face_norm=gp_Vec(face_norm)
        
        display, start_display, add_menu, add_function_to_menu = init_display()
        
        display.DisplayShape(exposed_building,color=gray,transparency=0.5)
        
        display.DisplayShape(self._face,color='BLUE',transparency=0.5)
        
        # face norm
        display.DisplayVector(face_norm,center)
        
        # Set the position to plot the vector representing
        # the direction of discretization
        lvec=[gp_Vec(d) for d in self._lvec]
        pos =[ center.Translated(v*(-3.)) for v in lvec]
        [display.DisplayVector(v,p) for v,p in zip(lvec,pos)]
        
        
        display.FitAll()
        start_display()

      
        


class direct_mask_on_face:
    def __init__(self,face,lsun_dir):
        self._face=face
        self._lsun_dir=lsun_dir
        self._mask_faces=[]
    
    def compute_mask(self,exposed_building):
        """Compute the shaded area on a given face for 
        multiple solar position

        Parameters
        ----------
        exposed_building : shape
        min_are : float
        Returns
        -------
        None
            
        """        
        srf = BRep_Tool().Surface(self._face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
        
        for j,sun_dir in enumerate(self._lsun_dir):
           
            mask_face=compute_direct_mask_on_face(sun_dir,exposed_building,self._face,face_norm)
            
            print('\r     sun dir ',j,'/',len(self._lsun_dir)-1,end="",flush=True)
            self._mask_faces.append(mask_face)
        
    def display(self,exposed_building):
        """Display shapes

        Parameters
        ----------
        exposed_building : Shape
        
        Returns
        -------
        None    
        """        
        
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)

        srf = BRep_Tool().Surface(self._face)
        umin,umax,vmin,vmax=breptools_UVBounds(self._face)
        center = srf.Value((umax+umin)*.5,(vmax+vmin)*.5)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
        face_norm=gp_Vec(face_norm)
        
        lvec=[gp_Vec(d) for d in self._lsun_dir]
        # display of sun vectors around the face
        sun_pos =[ center.Translated(v*(-3.)) for v in lvec]
        
        
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.DisplayVector(face_norm,center)
        [display.DisplayVector(v,p) for v,p in zip(lvec,sun_pos)]
        display.DisplayShape(exposed_building,color=gray,transparency=0.5)
        display.DisplayShape(self._face,color='BLUE',transparency=0.5)
        
        display.FitAll()
        start_display()
        
        
    
    
    def shadow_areas(self):
        gpp=GProp_GProps()
        l=[]
        for f in self._mask_faces:
            brepgprop_SurfaceProperties(f,gpp)
            l.append(gpp.Mass())
        
        self._mask_area=np.array(l)
        
    def face_area(self):
        gpp=GProp_GProps()
        brepgprop_SurfaceProperties(self._face,gpp)
        self._face_area=gpp.Mass()
        
    def mask_ratio(self):
        self._mask_ratio= 1.-self._mask_area/self._face_area
        
    def compute_complementary_face(self):
        """ Compute the shape of lighted part of a face

        Parameters
        ----------
        None
        -------
        The cutted building : Solid
        None    
        """        

        cutter=BOPAlgo_BOP()
        for f in self._mask_faces:
            if not f.IsNull():
                cutter.Clear()
                cutter.SetOperation(BOPAlgo_Operation.BOPAlgo_CUT)
                cutter.AddArgument(self._face)
                cutter.AddTool(f)
                cutter.Perform()
                self._complementary=cutter.Shape()
                #print(' cutter ',complementary)
                        
            else :
                
                self._complementary=complementary=glass_face  


class rtaa_on_faces:
    """Compute the RTAA shading factor (CM) on a given set of faces.
            
    """
    def __init__(self,lfaces,wall_plane,exposed_building,l_sun_dir,cut=True,adjust=True):
        
        self._lfaces=lfaces
        self._wall_plane=wall_plane
        #self._original_lfaces=[f for f in lfaces]
        self._exposed_building=exposed_building
        self._lsun_dir=l_sun_dir
        self._cut = cut
        self._adjust = adjust
        
    
    
    
    def cut_exposed_building(self,cut_size,input_face):
        """Reduce the size of building solid by cutting it.
        The face of interest (face_ref) is enlarged (cut_size) 
        and extruded trough the exterior normal. The building parts 
        intersecting the created volume will be used to 
        compute shadow

        Parameters
        ----------
        cut_size : float
        face_ref : Face
        Returns
        -------
        The cutted building : Solid
            
        """        
        
        
        srf = BRep_Tool().Surface(input_face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(input_face.Orientation()==1):
            face_norm.Reverse()
        
        pln=plane.Pln().Translated( gp_Vec(face_norm)*(-1))    
        newface= BRepBuilderAPI_MakeFace(plane.Pln(),-cut_size,cut_size,-cut_size,cut_size).Face()
        extrusion = BRepPrimAPI_MakePrism(newface,gp_Vec(face_norm)*10.,False,False).Shape()
        
        intersector=BOPAlgo_BOP()
        intersector.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
        intersector.AddTool(extrusion) 
        intersector.AddArgument(self._exposed_building)
        intersector.Perform()
        return shapes_as_solids([intersector.Shape()])[0]
    
    def adjust_face_to_wall_plane(self,input_face):
        """Translate the input face to be coplanar with the external 
        face of the hosting wall

        Parameters
        ----------
        input_face : Face
        Returns
        -------
        Translated face : Face
            
        """      
        ext_plane=self._wall_plane.Pln()
        srf = BRep_Tool().Surface(input_face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(input_face.Orientation()==1):
            face_norm.Reverse()
        
        plane=plane.Pln()
        plane_diff=gp_Vec(ext_plane.Location().XYZ())-gp_Vec(plane.Location().XYZ())
        distance = plane_diff.Dot(gp_Vec(face_norm))
        
        
        translation = gp_Trsf()
        translation.SetTranslation(gp_Vec(face_norm)*distance)
        builder = BRepBuilderAPI_Transform(translation)
        builder.Perform(input_face)
        #print("Adjsuted face : new face with translation ",distance)
       
        return builder.ModifiedShape(input_face)
        
    
    def compute_masks(self):
        """Cut the building, adjust the face positions, 
        compute direct and diffuse masks

        Parameters
        ----------
        None
        Returns
        -------
        None
            
        """        
        self._ldirect=[]
        self._ldiffuse=[]
        
        self._lcuted=[]
        self._ladjusted=[]
        
        for f in self._lfaces:
         #rtaa_on_faces(exposed_building,sun_to_earth_project)
            if self._cut :
                cutted_building = self.cut_exposed_building(1.,f)
            else :
                cutted_building = self._exposed_building
             
            if self._adjust :
                adjusted_face = self.adjust_face_to_wall_plane(f)
            else:
                adjusted_face = f
                
            print('direct mask')
            direct_mask = direct_mask_on_face(adjusted_face,self._lsun_dir)
            direct_mask.compute_mask(cutted_building)
            direct_mask.shadow_areas()
            direct_mask.face_area()
            direct_mask.mask_ratio()
            #direct_mask.display(cutted_building)
            print('\n')
            
            print('diffuse mask ')
            diffuse_mask=diffuse_mask_on_face(adjusted_face)
            diffuse_mask.compute_mask(cutted_building)
            diffuse_mask.compute_weighted_mask()
            print('\n')
            
            self._ldirect.append(direct_mask)
            self._ldiffuse.append(diffuse_mask)
            
            self._lcuted.append(cutted_building)
            self._ladjusted.append(adjusted_face)
    
    def display(self):
        """Display the part of the (eventually cutted) building used to compute 
        the shading factor

        Parameters
        ----------
        None
        Returns
        -------
        None
            
        """        
        
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)
    
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.DisplayShape(self._exposed_building,color='GREEN',transparency=0.5)
        
        [ display.DisplayShape(b,color=gray,transparency=0.5) for b in self._lcuted]
        gpp=GProp_GProps()
        for f in self._ladjusted:
            display.DisplayShape(f,color='RED',transparency=0.5)
            
            srf = BRep_Tool().Surface(f)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(f.Orientation()==1):
                face_norm.Reverse()
            norm=gp_Vec(face_norm)
            brepgprop_SurfaceProperties(f,gpp)
            mc=gpp.CentreOfMass()
            display.DisplayVector(norm,mc)
        
        
        display.FitAll()
        
        start_display()
     
    
    
    def compute_masks_hangover(self,hp,lp,dhp,dhm,pcd=0.0,pcg=0.0):
        """Reproduce partially the value of CM computed by excel sheet
           Not complete !!

        Parameters
        ----------
        
        Returns
        -------
        None
            
        """
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        Fdir=np.array(len(self._lsun_dir))
        for face in self._lfaces:
            srf = BRep_Tool().Surface(face)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(face.Orientation()==1):
                face_norm.Reverse()
            
                            
            cosbeta= Zaxis.Direction().Dot(face_norm)
            Zvec = Zaxis.Direction()
            earth_to_sun =  [s.Reversed() for s in self._lsun_dir]
            zenith  = np.array([Zaxis.Direction().Angle(s) for s in earth_to_sun])
            directflux = np.array([face_norm.Dot(s) for s in earth_to_sun])
            horizontalflux= np.array([Zvec.Dot(s) for s in earth_to_sun])
            
            Zvec=gp_Vec(Zaxis.Direction())
            svec= [ gp_Vec(s) for s in earth_to_sun]
            projXY = [Zvec.CrossCrossed(s,Zvec) for s in svec]
            gamma = [xy.Angle(s) for s,xy in zip(svec,projXY)]
            psialpha=[xy.Angle(gp_Vec(face_norm)) for  xy in projXY]
            
            dhmprime=dhm*np.tan(gamma)/np.cos(psialpha)
            gammaprime=np.sin(psialpha)/np.tan(gamma)
            lfdir=[]
            for gammap,dhmp in zip(gammaprime,dhmprime):
                A1d=np.max([ 0, np.min([lp ,  lp+pcd    +np.min([dhp,dhmp])*gammap ])])
                A2d=np.max([ 0, np.min([hp ,  dhmp-dhp, -pcd/gammap -dhp])])
                A3d=np.max([ 0, np.min([lp ,  lp+pcd    +np.min([dhp+hp,dhmp])*gammap ])])
                A4d=np.max([ 0, np.min([hp ,  dhmp-dhp, -(pcd+lp)/gammap -dhp ])])
                Fdird=1-(A1d*A2d + (A4d-A2d)*(A1d+A3d)*.5)/(hp*lp)
                
                A1g=np.max([ 0, np.min([lp ,  lp+pcg   -np.min([dhp,dhmp])*gammap ])])
                A2g=np.max([ 0, np.min([hp ,  dhmp-dhp, pcg/gammap -dhp])])
                A3g=np.max([ 0, np.min([lp ,  lp+pcg   -np.min([dhp+hp,dhmp])*gammap ])])
                A4g=np.max([ 0, np.min([hp ,  dhmp-dhp, (pcg+lp)/gammap -dhp ])])
                Fdirg=1-(A1g*A2g + (A4g-A2g)*(A1g+A3g)*.5)/(hp*lp)
            
                lfdir.append(Fdirg*Fdird)
            
            not_exposed= directflux<1e-5
            
            self._LFdir=np.array(lfdir)
            self._LFdir[not_exposed]=0.0
            
            
            alpha= np.arctan((dhp+hp*.5)/dhm)
            beta = np.arctan((dhp+.5*hp)/(lp*.5+(pcg+pcd)*.5))
            self._Fdiff=alpha/(.5*np.pi)*(1.- beta/(np.pi/2.))+beta/(.5*np.pi)
        
    def compute_cm(self,data_irradiance):
        """Compute the shading factor according to the RTAADOM

        Parameters
        ----------
        data_irradiance : DataFrame
            Dataframe of DHI, DNI and GHI at each solar position
        
        Returns
        -------
        None
            
        """
        albedo=.2
        origin = gp_Pnt(0.0,0.0,0.0)
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        
        self._cm_byface=[]
           
        self._ldf_irr=[]
        
        for face,dir,diff in zip(self._lfaces,self._ldirect,self._ldiffuse):
            
            srf = BRep_Tool().Surface(face)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(face.Orientation()==1):
                face_norm.Reverse()
            
                            
            cosbeta= Zaxis.Direction().Dot(face_norm)
            Zvec = Zaxis.Direction()
            earth_to_sun =  [s.Reversed() for s in self._lsun_dir]
            zenith  = np.array([Zaxis.Direction().Angle(s) for s in earth_to_sun])
            directflux = np.array([face_norm.Dot(s) for s in earth_to_sun])
            horizontalflux= np.array([Zvec.Dot(s) for s in earth_to_sun])
            
            not_exposed = directflux<-1e-5
            directflux [not_exposed]=0.0
            
            Rb = directflux/horizontalflux
            
            maskratio = dir._mask_ratio
                        
            Drp = data_irradiance['DNI']*directflux
            Dfp = data_irradiance['DHI']*(1+cosbeta)*.5
            Rrp = data_irradiance['GHI']*albedo*(1-cosbeta)*.5
            
            # RTAA rule
            Drp[not_exposed]=0.0
            #Rrp[not_exposed]=0.0
            
            
            
            mDrp = Drp*maskratio
            mDfp = Dfp*diff._wm
            mRrp = Rrp*diff._wm_soil
            
            masked_irr   = mDrp + mDfp + mRrp 
            unmasked_irr =  Drp +  Dfp +  Rrp
            
            self._cm_dir = mDrp.sum()/Drp.sum()
            self._cm_diff= diff._wm
            self._cm_ref = diff._wm_soil
            
            self._cm_byface.append( masked_irr.sum()/unmasked_irr.sum() )
                        
            face_irr=dict()
            face_irr['DNI']=data_irradiance['DNI']
            face_irr['DHI']=data_irradiance['DHI']
            face_irr['GHI']=data_irradiance['GHI']
            
            face_irr['Drp']=Drp
            face_irr['Dfp']=Dfp
            face_irr['Rrp']=Rrp
            
            face_irr['mDrp']=mDrp
            face_irr['mDfp']=mDfp
            face_irr['mRrp']=mRrp
            
            face_irr['mask']=maskratio
            face_irr['angle']=np.array([face_norm.Angle(s) for s in earth_to_sun])
            face_irr['cm_inst']=masked_irr/unmasked_irr
            
            self._ldf_irr.append(pd.DataFrame(face_irr))
           
                   
        
        all_values=pd.concat(self._ldf_irr,axis=0)
        total_irr =all_values[['Drp','Dfp','Rrp']].sum().sum()
        total_m_irr =all_values[['mDrp','mDfp','mRrp']].sum().sum()
        
        self._cm= total_m_irr/total_irr


    
class rtaa_solar_study:
    """The object that manage the data to compute shading factor on
    selected faces of a building described by a ifc file
        
    """
    def __init__(self,ifcfilename):
        
        self._ifc_file= ifcopenshell.open(ifcfilename)
        self._solar_elements=dict()
        self._building_elements=dict()
        print(" ifc file ",self._ifc_file)
        
        self._proj_loc=project_location()
        self._proj_loc.set_location_from_ifc(self._ifc_file)
        self._proj_loc.set_northing_from_ifc(self._ifc_file)
        self._proj_loc.load_irradiance()
    
        
    def _add_elements(self,ids=[],ltypes=[],container=None):
        if((len(ids)==0) & (len(ltypes)==0)):
            print(" Error : provides ids or classes to set elements of interest")
            return
        if container is None:
            print(" Error : provide a container")
            
        element_set=set()
        for id in ids:
            element_set.add(self._ifc_file.by_id(id))
        
        for t in ltypes:
            [element_set.add(el) for el in self._ifc_file.by_type(t)]
          
        
        if container=="solar":
            self._solar_elements={el.id():el for el in element_set}
            
        elif container=="building":
            self._building_elements={el.id():el for el in element_set}
        
    
    
    def _remove_elements_by_ids(self,ids=[],container=None):
        if(len(ids)==0):
            print(" Error : provides ids or classes to set elements of interest")
            return
        if container is None:
            print(" Error : provide a container")
            
                    
        if container=="solar":
            for id in ids:
                self._solar_elements.pop(id)
        elif container=="building":
            for id in ids:
                self._building_elements.pop(id)
    
    
    def add_solar_elements(self,ids=[],types=[]):
        self._add_elements(ids,types,'solar')
        print(" Number of elements for solar analysis: ",len(self._solar_elements.keys()))
    
    def remove_solar_elements(self,ids=[]):
        self._remove_elements(ids,types,'solar')
        
    def add_building_elements(self,ids=[],types=[]):
        self._add_elements(ids,types,'building')
        print(" Number of elements in building: ",len(self._building_elements.keys()))
    
    def remove_building_elements(self,ids=[]):
        self._remove_elements(ids,types,'building')    
                
        
    def set_geometries(self):
        self._building_solid()
        self._solar_faces()
                
    
    def _solar_faces(self):
        """

        Parameters
        ----------
        
        Returns
        -------
        None
            
        """
        self._solar_shapes = { x.id():ifcelement_as_solid(x) 
                            for x in self._solar_elements.values() }
        
        self._solar_faces=defaultdict(list)
        self._hosting_solar_id=dict()
        
        for id,shape in self._solar_shapes.items():
             
            if self._solar_elements[id].is_a('IfcWall'):
                hosting_wall_id=id
            elif((self._solar_elements[id].is_a('IfcWindow')) | 
                  (self._solar_elements[id].is_a('IfcDoors'))):
                el = self._solar_elements[id]
                citedby=list(self._ifc_file.get_inverse(el.FillsVoids[0].RelatingOpeningElement))
                for c in citedby:
                    if c.is_a('IfcRelVoidsElement'):
                        hosting_wall=c.RelatingBuildingElement
                        hosting_wall_id = hosting_wall.id()
            else:
                print("Error, solar element is not a wall,a  window or a door")
            self._hosting_solar_id[id]=hosting_wall_id
            
            if(self._wall_norm_plane[hosting_wall_id]):
                wall_norm,plane = self._wall_norm_plane[hosting_wall_id]
                
                faces = biggestface_along_vector(shape,wall_norm)
                if len(faces)>0:
                    self._solar_faces[id]=faces
            else :
                print(" No normal found for this window ",id)
        

        
    def _building_solid(self):
        """Set the geometry of building according to the selected elements

        Parameters
        ----------
        None

        Returns
        -------
        None
            
        """
        building_shapes = { x.id():ifcelement_as_solid(x) 
                            for x in self._building_elements.values() }
        
        self._building= fuse_listOfShape(list(building_shapes.values()))
        
        # # Compute the external shell of the building
        # # Fill the wall building with the spaces as solid
        
        # get spaces geometries
        ifcspaces=self._ifc_file.by_type('IfcSpace')
        # list of boolean
        in_building=[ space.id() not in list(self._building_elements.keys()) for space in ifcspaces]
        # space in building
        space_shapes = [ifcelement_as_solid(x) for x in compress(ifcspaces,in_building) ]
        # making list of solids
        solids=shapes_as_solids([self._building]+space_shapes)
        
        # compute the external shell
        self._external_shell = get_external_shell(solids)
        
        
        wall_shapes={ el.id():building_shapes[el.id()] for el in self._building_elements.values() if el.is_a()=='IfcWall'}
        self._wall_norm_plane={}
        for id,wall_shape in wall_shapes.items():
            self._wall_norm_plane[id]=exterior_wall_normal_and_plane(wall_shape,self._external_shell)
        
    def display(self,disp_external_shell=False,disp_solar_shapes=False):
        """Display the building, shading factor faces and external enveloppe

        Parameters
        ----------
        disp_external_shell : bool, optional
            Display the computed external envelop of the building
        disp_solar_shapes : bool, optional
            Display the faces on which shading factor is computed

        Returns
        -------
        None
            
        """
        
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)
    
        display, start_display, add_menu, add_function_to_menu = init_display()
    
        display.DisplayShape(self._building,color=gray,transparency=0.5)
        
        if disp_external_shell:
            display.DisplayShape(self._external_shell,color='BLUE',transparency=0.5)
            
        if disp_solar_shapes:
            [display.DisplayShape(s,color='GREEN',transparency=0.5) for s in self._solar_shapes.values()]
        
        gpp=GProp_GProps()
        
        for id,facelist in self._solar_faces.items():
            idhost = self._hosting_solar_id[id] 
            norm,_ =self._wall_norm_plane[idhost]
            norm=gp_Vec(norm)
            for f in facelist:
                display.DisplayShape(f,color='RED',transparency=0.1)
                
                brepgprop_SurfaceProperties(f,gpp)
                mc=gpp.CentreOfMass()
                display.DisplayVector(norm,mc)
        
        display.FitAll()
        
        start_display()
        
        
    
    def run(self,cut=True,adjust=True):
        """Make the computation of shading factor for the selected openings

        Parameters
        ----------
        cut : bool, optional
            Cut the building to operate with smaller model
        adjust : bool, optional
            Make glass faces coplanar with hosting wall external wall face 

        Returns
        -------
        None
            
        """

        
        self._results=dict()
        
        for id,facelist in self._solar_faces.items():
            
            
            idhost = self._hosting_solar_id[id] 
            norm,plane =self._wall_norm_plane[idhost]
            
            sun_to_earth_project,irradiance=self._proj_loc.data_critical_period_day(facelist[0])
            
            print(" Opening ID  : ", id)
            rtaa=rtaa_on_faces(facelist,plane,self._building,sun_to_earth_project,cut,adjust)
            
            rtaa.compute_masks()
            rtaa.compute_cm(irradiance)
            self._results[id]=rtaa
            #rtaa.display()
            
        
    def cm(self):
        return { id : r._cm for id,r in self._results.items()} 
        
    
        
if __name__ == "__main__":

    
    filename = '../tests/data/debords_casquettes_fins.ifc'
    
    rss=rtaa_solar_study(filename)
    rss.add_building_elements([],['IfcWall','IfcSlab'])
    rss.add_solar_elements([],['IfcWindow'])
    rss.set_geometries()
    rss.run()