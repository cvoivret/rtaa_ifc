from collections import defaultdict
import itertools
from array import array
import numpy as np
import pandas as pd
import ifcopenshell
import sunposition as sunpos
from  datetime import datetime
from timezonefinder import TimezoneFinder
import pytz

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.tri as tri

import gc
import scipy.integrate as spi

import line_profiler

from timeit import default_timer as timer

from ifcopenshell.geom import create_shape
from ifcopenshell.geom.occ_utils import yield_subshapes

from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color,Quantity_TOC_RGB

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox,BRepPrimAPI_MakePrism,BRepPrimAPI_MakeHalfSpace,BRepPrimAPI_MakeSphere,BRepPrimAPI_MakeCylinder
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties,brepgprop_VolumeProperties,brepgprop_LinearProperties
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing,BRepBuilderAPI_MakeSolid,BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Copy,	BRepBuilderAPI_Transform
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepTools import breptools_UVBounds

from OCC.Core.GProp import GProp_GProps
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.gp import gp_Pnt,gp_Dir,gp_Vec,gp_Pln,gp_Lin,gp_Trsf,gp_Ax3,gp_Ax1
from OCC.Core.Geom import Geom_Plane

from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.TopTools import TopTools_ListOfShape,TopTools_IndexedMapOfShape
from OCC.Core.TopExp import topexp_MapShapes
from OCC.Core.TopAbs import TopAbs_SOLID,TopAbs_FACE,TopAbs_SHELL,TopAbs_WIRE

from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
 
from OCC.Core.BOPAlgo import BOPAlgo_BOP,BOPAlgo_Operation
from OCC.Core.BOPAlgo import BOPAlgo_CellsBuilder
from OCC.Core.BOPTools import BOPTools_AlgoTools_OrientFacesOnShell

from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib

import OCC.Core.ShapeFix as ShapeFix_Shape

from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector

from OCC.Core.Standard import standard_Purge

from OCC.Extend.DataExchange import write_stl_file


def fuse_listOfShape(los,FuzzyValue=1e-6):
    """
    Simple funtion to wrap boilerplate code
    """    
    fuser=BOPAlgo_BOP()
    fuser.SetOperation(BOPAlgo_Operation.BOPAlgo_FUSE)
    los_1=TopTools_ListOfShape()
    [los_1.Append(s) for s in los[::2]]
    los_2=TopTools_ListOfShape()
    [los_2.Append(s) for s in los[1::2]]
    fuser.SetArguments(los_1)
    fuser.SetTools(los_2)
    #fuser.SetFuzzyValue(FuzzyValue)
    fuser.SetNonDestructive(True)
    fuser.Perform()
    return fuser.Shape()

def shapes_as_solids(lshape):
    """
    Try to build a list of solid from a list to shapes.
    Flatten nested shapes.
    Sew tesselated shape to build a solid.
    Fix shape if needed.
    
    """
    lsolid=[]
     
    maps=TopTools_IndexedMapOfShape()
    for s in lshape:
        maps.Clear()
        topexp_MapShapes(s,TopAbs_SOLID,maps)
        if(maps.Size()>0):
            lsolid.extend([maps.FindKey(i) for i in range(1,maps.Size()+1)])
        else:
            maps.Clear()
            topexp_MapShapes(s,TopAbs_FACE,maps)
            sewer=BRepBuilderAPI_Sewing()
            [sewer.Add(maps.FindKey(i)) for i in range(1,maps.Size()+1)]
            sewer.Perform()
            sewed=sewer.SewedShape()
            if(sewed.ShapeType()==0):
                lshell=list(yield_subshapes(sewed))
                
                for shell in lshell:
                    lsolid.append(BRepBuilderAPI_MakeSolid(shell).Solid())
            else:
                solid=BRepBuilderAPI_MakeSolid(sewed).Solid()
                lsolid.append(solid)
    lsolid2=[]            
    for s in lsolid:
        fixer=ShapeFix_Shape.ShapeFix_Shape(s)
        fixer.Perform()
        lsolid2.append(fixer.Shape())
         
    return lsolid2

def get_external_shell(lshape):
    """
    try to identigy a shell (set of face) that limit the inside and the outside of the building.
    Basically, wall and room must be part of the solid list input.
    Build the boundign box of the whole model and enlarge it.
    
    """
       
    #unionize solids
    unionsolid=fuse_listOfShape(lshape)
    
    # Create the bounding box of the unioned model
    box=Bnd_Box()
    brepbndlib.Add(unionsolid,box)
    box.Enlarge(1.) # to avoid parallel face of bbox with face model to be coplanar(
    # broke shell intersection 
    boxshape=BRepPrimAPI_MakeBox(box.CornerMin(),box.CornerMax()).Shape()

    #boolean difference between the unioned model and its bounding box
    diff=BOPAlgo_BOP()
    diff.SetOperation(BOPAlgo_Operation.BOPAlgo_CUT)
    diff.AddArgument(boxshape)
    diff.AddTool(unionsolid)
    #diff.SetFuzzyValue(1e-5)
    diff.Perform()
    diffshape=diff.Shape()
    

    # boolean common of shells : could be considered as the shell 
    # separating interior and exterior of the building
    top=TopologyExplorer(unionsolid)
    unionshell = top.shells()
    

    top=TopologyExplorer(diffshape)
    diffshell = top.shells()


    common=BOPAlgo_BOP()
    common.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
    args=TopTools_ListOfShape()
    for shell in unionshell:
        args.Append(shell)
        
    tools=TopTools_ListOfShape()
    for shell in diffshell:
        tools.Append(shell)
        
    common.SetArguments(args)
    common.SetTools(tools)
    #common.SetFuzzyValue(1e-5)
    common.Perform()
    commonshell=common.Shape()
    
    BOPTools_AlgoTools_OrientFacesOnShell(commonshell)
    
    return commonshell

#@profile

def shadow_caster_ext_on_face(sun_dir,mybuilding,myface,theface_norm,min_area = 1e-3):
    """
    sun_dir = one vector (downward direction)
    building = a solids that possibily make shadow on face
    face = a face to cast shadow on from building along sun_dir
    face_norm = pointing to the exterior of the face (outside)
    
    return  : a face with zero or positive area, None if no shadow
    
    """
    # direction pointing downward (nigth time)
    #print(" sun vec   ", sun_dir.Coord())
    #print(" face norm ",theface_norm.Coord())
    #print(" Dot       ",theface_norm.Dot(sun_dir))
    #print(sun_dir.Z())
    
    #lmod.append(building.Modified())
    
    building = BRepBuilderAPI_Copy(mybuilding).Shape()
    theface = BRepBuilderAPI_Copy(myface).Shape()

    
    #if sun_dir.Z()>0.:
        #print('Z positive')
    #    return theface
    
    
    # face not exposed to the sun
    #print(' DOT ' , theface_norm.Dot(sun_dir))
    if theface_norm.Dot(sun_dir)>-1e-5:
        #print('not exposed',flush=True)
        #print('Not oriented to the sun')
        return theface
    
    gpp=GProp_GProps()
    brepgprop_SurfaceProperties(theface,gpp)
    gf_area=gpp.Mass()
    
    ext_vec=gp_Vec(sun_dir)
    ext_vec.Multiply(10.)
    
    # extrusion 
    extrusion1=BRepPrimAPI_MakePrism(theface,-ext_vec,False,True).Shape()
    #print("type ", extrusion1.
    """
    common=BOPAlgo_BOP()
    common.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
    common.AddArgument(extrusion1)
    common.AddTool(building)
    common.Perform()
    intersection=common.Shape()
    warnings=cb.DumpWarningsToString()
    warnings=warnings.splitlines()
    print(warnings)
    """
    
    los2 = TopTools_ListOfShape()
    los_avoid = TopTools_ListOfShape()
    los2.Append(building)
    los2.Append(extrusion1)
    
    #warnings = True
    precision = 1e-9
    ntry=1
    haswarning=True
    
    cb=BOPAlgo_CellsBuilder()
    cb.SetNonDestructive(True)
    cb.SetArguments(los2)
    cb.SetFuzzyValue(precision)

    #cb.SetGlue(True)
    cb.Perform()
    cb.RemoveAllFromResult()
    cb.AddToResult(los2,los_avoid,2)
    
    cb.RemoveInternalBoundaries()
    haswarning =   cb.HasWarnings()
    
    while( haswarning ):
        precision*=1.2
        ntry+=1
        #print(" precision ", precision)
                    
        cb=BOPAlgo_CellsBuilder()
        cb.SetNonDestructive(True)
        cb.SetArguments(los2)
        cb.SetFuzzyValue(precision)

        #cb.SetGlue(True)
        cb.Perform()
        cb.RemoveAllFromResult()
        cb.AddToResult(los2,los_avoid,2)
        
        cb.RemoveInternalBoundaries()
        haswarning =   cb.HasWarnings()
        #print(haswarning)
        #print(cb.DumpWarningsToString())
        
        if ( ntry>50):
            break
            
        
        
        #print(" precision update ", precision)
        
    #lntry.append(ntry)    
    #lprec.append(precision)
        
    intersection=cb.Shape()
    #h=cb.History()
    
        
    #print('hist ',h.DumpToString())
    
    """
    common=BRepAlgoAPI_Common(building,extrusion1)
    extrusion1=common.Shape()
    """
    #print("first")
    
    """
    if len(warnings)>0 and len(warnings) <20 :
    
        print("Warning ",warnings)
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)

        display, start_display, add_menu, add_function_to_menu = init_display()

        display.DisplayShape(building,color=gray,transparency=0.9)
        display.DisplayShape(extrusion1,color='RED',transparency=0.9)

        display.DisplayShape(intersection,color='BLUE',transparency=0.1)
    
        display.FitAll()
        #ifcopenshell.geom.utils.main_loop()
        start_display()
    """
    
    #print("Error ",common.DumpErrorsToString())
 
 
    intersection_faces=list(TopologyExplorer(intersection).faces())
    if(len(intersection_faces)==0):
        
        return TopoDS_Face()           
    # flux trough the face
    brepgprop_SurfaceProperties(theface,gpp)
    theface_area=gpp.Mass()
    theface_flux=theface_norm.Dot(sun_dir)*theface_area
    #print("face area ",theface_area)
    #print("face_flux ",theface_flux)
    
    larea=[]
    lfaces=[]
          
    for ff in intersection_faces:
        
        adapt=BRepAdaptor_Surface(ff)
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
            if(ff.Orientation()==1):
                ff.Reverse()
            faceflux = larea[-1]*fn.Dot(sun_dir)
            #print("   unitflux ",faceflux)
            if (faceflux/theface_flux>1e-4):
                lfaces.append(ff)
            
            #lfaces.append(ff)
    
    lsolid=[ BRepPrimAPI_MakePrism(s,ext_vec,False,True).Shape() for s in lfaces]
    
    
    if(len(lsolid)==0):
        return TopoDS_Face() # void face with zero area
    
    brepgprop_SurfaceProperties(theface,gpp)
    totarea=gpp.Mass()
    
    lface2=[]
    for s,f in zip(lsolid,lfaces):
        common=BRepAlgoAPI_Common(s,theface)
        sh=common.Shape()
        """
        brepgprop_SurfaceProperties(sh,gpp)
        area_proj=gpp.Mass()
        #brepgprop_SurfaceProperties(f,gpp)
        #area=gpp.Mass()
        if(area_proj/totarea<1e-4):
            continue
        """
        lface2.append(sh)
    #print(" total face ",len(lfaces))    
    #print(" with intersection :",len(lface2)    )
    #print(" with flux :",len(lfaces3))
    
    if len(lface2)==1:
        shadowface=lface2[0]
    else:    
        los2 = TopTools_ListOfShape()
        for s in lface2:
            los2.Append(s)
        
        cb=BOPAlgo_CellsBuilder()
        cb.SetArguments(los2)
        cb.Perform()
        cb.AddAllToResult(2,False)
        cb.RemoveInternalBoundaries()
        shadowface=cb.Shape()
    
    if shadowface==None:
        shadowface=TopoDS_Face()
    
    
    return shadowface


def compute_diffuse_mask_on_face(mybuilding,theface,theface_norm,min_area):
    """
    Compute the matrix of mask for isotropic diffuse shading
    must be computed only once for a face
    """
    
    # set a discretization of the half sphere exposed to the sun
    n=10
    
    arct=np.arccos(np.linspace(0,1,n)[::-1])
    
    temp=np.arccos(np.linspace(1,-1,n))
    arcp=np.concatenate((temp,temp[1:]+np.pi))
    
    #t= np.linspace(0,.5*np.pi,n  ,endpoint=True)
    #p= np.linspace(0, 2*np.pi,2*n,endpoint=True)
    t=arct
    p=arcp
    """
    print('arct ',arct)
    print(' t',t)
    print('\narcp ',arcp)
    print(' p',p)
    cdsc
    """
    mid_t=.5*(t[1:]+t[:-1])
    mid_p=.5*(p[1:]+p[:-1])
    # coordinate to set the direction of mask calculation
    t1,p1=np.meshgrid(mid_t,mid_p,indexing='ij')
    #surf = np.zeros((t.shape[0]-1,p.shape[0]-1))
    areas = np.zeros_like(t1)
    
    # computation of local area associated with each direction
    # sum = 2*pi
    
    ds=lambda x,y: np.sin(x)
    
    for i,j in itertools.product(range(t1.shape[0]),range(t1.shape[1])):    
        
        res,err = spi.nquad(ds, [[t[i], t[i+1]],[p[j], p[j+1]]])
        areas[i,j]=res
    areas=areas.flatten()   
    
    # first incline (between 0+delta (normal) and PI/2) with respect of one of the axis of the surface 
    # we dont care which axis because of the next rotation
    srf = BRep_Tool().Surface(theface)
    adapt=BRepAdaptor_Surface(theface)
    #umin,umax,vmin,vmax=breptools_UVBounds(ff)
    #props=GeomLProp_SLProps(srf3,0.5*(umax-umin),0.5*(vmax-vmin),1,0.001)
    #fn=props.Normal()
    plane = adapt.Plane()
    Xaxis_plane_ax=plane.XAxis()
    normal_plane_ax = plane.Axis()
    normal_plane_vec = gp_Vec(normal_plane_ax.Direction()) # unit vector
    
    origin = gp_Pnt(0.0,0.0,0.0)
    Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
    
    #print('Xaxis plane ',Xaxis.Direction().Coord())
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
        lmask_sky.append( final.Dot(gp_Vec(Zaxis.Direction()))>0.)
        #lvector_diff.append(final)
        
        projected = normal_plane_vec.CrossCrossed(final,normal_plane_vec)
        #print(projected.Magnitude())
        lvec_proj.append(projected.Coord())
    
    #print(' normal ',theface_norm.Coord())
    #print(' transformed ', transformed.Coord())
    #print(lmask_sky)
    
    # remove the axis with null coordinate (in surface plan)
    projected_vec= np.array(lvec_proj)
    #print(projected_vec)
    idx=np.argwhere(np.all(np.isclose(projected_vec[...,:],np.zeros(projected_vec.shape)),axis=0))
    #print(idx)
    projected_vec=np.delete(projected_vec,idx,axis=1)
    
        
    # rotate with respect to the face normal 
    
    # negate the vectors
    
    
    
    # compute ratio for each direction
    dmof=direct_mask_on_faces([theface],lvec)
    dmof.compute_mask(mybuilding,min_area)
    dmof.compute_area_and_ratio()
    
    mask_vector = np.array(dmof._ratio_vector)
    mask_sky = np.array(lmask_sky)
    
    return t1,p1,areas,mask_vector,projected_vec,mask_sky 
    
    # could use sof with discretization vector instead of sun_dir
    # Sum up and return


 
#@profile
def shadow_caster_ext2(sun_dir,building,theface,theface_norm,min_area = 1e-3):
    """
    sun_dir = one vector (downward direction)
    building = a solids that possibily make shadow on face
    face = a face to cast shadow on from building along sun_dir
    face_norm = pointing to the exterior of the face (outside)
    
    return  : a face with zero or positive area, None if no shadow
    
    """
    #print(theface_norm.Dot(sun_dir))
    # face not exposed to the sun
    if theface_norm.Dot(sun_dir)>-1e-5:
        #print('not exposed',flush=True)
        return theface# void face with zero area
    gpp=GProp_GProps()
    brepgprop_SurfaceProperties(theface,gpp)
    gf_area=gpp.Mass()
    
    ext_vec=gp_Vec(sun_dir)
    ext_vec.Multiply(2)
    
    # extrusion of 
    extrusion1=BRepPrimAPI_MakePrism(theface,-ext_vec,False,True).Shape()
    
    """
    intersector=BOPAlgo_BOP()
    intersector.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
    intersector.AddTool(extrusion1) 
    intersector.AddArgument(building)
    intersector.Perform()
    intersection=intersector.Shape()
    """
    common1=BRepAlgoAPI_Common(extrusion1,building)
    #common1.Build()
    intersection=common1.Shape()

    
    intersection_faces=list(TopologyExplorer(intersection).faces())
               
    larea=[]
    lfaces=[]
    
    for ff in intersection_faces:
        
        adapt=BRepAdaptor_Surface(ff)
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
            if(ff.Orientation()==1):
                ff.Reverse()
            
            lfaces.append(ff)
    
    lsolid=[ BRepPrimAPI_MakePrism(s,ext_vec,False,True).Shape() for s in lfaces]
    
    
    if(len(lsolid)==0):
        return TopoDS_Face() # void face with zero area
    
    brepgprop_SurfaceProperties(theface,gpp)
    totarea=gpp.Mass()
    
    """
    lface2=[]
    

    for s,f in zip(lsolid,lfaces):
        
    
        common=BRepAlgoAPI_Common(s,theface)
        
        sh=common.Shape()
        
        brepgprop_SurfaceProperties(sh,gpp)
        area_proj=gpp.Mass()
        #brepgprop_SurfaceProperties(f,gpp)
        #area=gpp.Mass()
        if(area_proj/totarea<1e-4):
            continue
        lface2.append(sh)
    """
    ## keep that, remove intersection with unit face, intersect with the union
    # profile !!
    los3 = TopTools_ListOfShape()   
    [los3.Append(s) for s in lsolid]    
    cb1=BOPAlgo_CellsBuilder()
    cb1.SetArguments(los3)
    cb1.Perform()
    cb1.AddAllToResult(2,False)
    cb1.RemoveInternalBoundaries()
    shape2=cb1.Shape()
    #lshapee.append(shape2)    
    
    if not shape2:
        return TopoDS_Face()
        #len(lface2)==1:
        #shadowface=lface2[0]
    else:    
        
        common=BRepAlgoAPI_Common(shape2,theface)
        
        shadowface=common.Shape()
        
    
    if shadowface==None:
        return TopoDS_Face()
       
    return shadowface





def shadow_caster_ray(sun_dir,building,theface,theface_norm,Nray=5):
 
    sphere_rad=0.05
    lshape=[]
        
    #discretize the face with Nray points
    srf = BRep_Tool().Surface(theface)
    umin,umax,vmin,vmax=breptools_UVBounds(theface)
    
    uoffset=0.5*(umax-umin)/Nray
    voffset=0.5*(vmax-vmin)/Nray
    
    uvalues,vvalues= np.meshgrid(np.linspace(umin+uoffset,umax-uoffset,Nray),
                                 np.linspace(vmin+voffset,vmax-voffset,Nray))
    
    # face not exposed to the sun
    if theface_norm.Dot(sun_dir)>-1.e-5:
        
        for u,v in zip(uvalues.flatten(),vvalues.flatten()):
            point=srf.Value(u,v)
            #lshape.append(BRepPrimAPI_MakeSphere(point,sphere_rad).Shape())
        return np.ones(uvalues.shape)#,lshape# all points of discretization are in shadow
    
    
    shape_inter = IntCurvesFace_ShapeIntersector()
    shape_inter.Load(building, 1e-6)
    infinity=float("+inf")
    nbpoints=array('b')
    for u,v in zip(uvalues.flatten(),vvalues.flatten()):
        point=srf.Value(u,v)
        line=gp_Lin(point,-sun_dir)
        shape_inter.PerformNearest(line, 0.0,100.)
        nbpoints.append(shape_inter.NbPnt())
        #if(shape_inter.NbPnt()>0):
        #    lshape.append(BRepPrimAPI_MakeSphere(point,sphere_rad).Shape())
    
    #print(nbpoints)
    res=np.array(nbpoints).reshape(uvalues.shape)
    
    res[res>0.]=1
    return res 

def exterior_wall_normal_and_plane(wall_shape,external_shell):
    
    gpp=GProp_GProps()
    
    tools=TopTools_ListOfShape()
    tools.Append(external_shell) 

    args=TopTools_ListOfShape()
    args.Append(wall_shape)
        
    common=BOPAlgo_BOP()
    common.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
        
    common.SetArguments(args)
    common.SetTools(tools)
    
    #common.SetFuzzyValue(1e-6)
    common.Perform()
    commonshell2=common.Shape() 
        
    
    # exteriro wall !!
    if commonshell2:
        faces=list(TopologyExplorer(commonshell2).faces())
        norm_area=defaultdict(float)
        norm_map=defaultdict(list)
        for f in faces:
            srf = BRep_Tool().Surface(f)
            plane = Geom_Plane.DownCast(srf)
            fn = plane.Axis().Direction()
            if(f.Orientation()==1):
                fn.Reverse()
            face_norm_coord=fn.Coord()
            # maybe necessary to round...
            face_norm_coord = tuple(round(c,10) for c in face_norm_coord)
            brepgprop_SurfaceProperties(f,gpp)
            norm_area[face_norm_coord]+=gpp.Mass()
            norm_map[face_norm_coord].append(f)
        wall_norm = max(norm_area, key=norm_area.get)   
                
        #print(norm_area)
        #print(wall_norm)
        #print(norm_map[wall_norm])
                    
        # wall_norm is rounded but almost equal to all element in the list
        # taking the first
        #lface_wall.append(norm_map[wall_norm])
        #lface_wall.append(norm_map[wall_norm][0])
        first_wall_face =norm_map[wall_norm][0]
        srf = BRep_Tool().Surface(first_wall_face)
        plane = Geom_Plane.DownCast(srf)
        wall_norm = plane.Axis().Direction()
        if(first_wall_face.Orientation()==1):
            wall_norm.Reverse()
        
        return wall_norm,plane
    

def exterior_wall_normal_dict(wallwindow,external_shell):
    
    wallnorm=dict()
    #shape of wall with a window
    wall_shapes=[create_shape(setting, ifc_file.by_guid(w_id)).geometry 
                    for w_id in wallwindow.keys() if ifc_file.by_guid(w_id).Representation 
                    is not None]
    for (w_id,ws) in zip(wallwindow.keys(),wall_shapes):
        wallnorm[w_id]=exterior_wall_normal(ws,external_shell)
    
    return wallnorm
    
def biggestface_along_vector(shape,vector,tol=1e-6,ratio=0.9):
    gpp=GProp_GProps()
    faces=list(TopologyExplorer(shape).faces())
    #print(" nb face par fenetre ", len(faceswin))
    facelist=[]
    facearea=[]
    #facenormal=[]
    for f in faces:
        top=TopologyExplorer(f)
        #print(top.number_of_wires())
        # face with some kind of hole
        if top.number_of_wires()>1:
            continue
        srf = BRep_Tool().Surface(f)
        plane2 = Geom_Plane.DownCast(srf)
        face_norm = plane2.Axis().Direction()
        if(f.Orientation()==1):
            face_norm.Reverse()
        
        
        if(face_norm.IsEqual(vector,tol)):
            #print(" face2 ",win_norm.Coord())
           
            brepgprop_SurfaceProperties(f,gpp)
            #print(" area ", gpp.Mass())
            facearea.append(round(gpp.Mass(),5))
            facelist.append(f)
            #facenormal.append(face_norm)
    #print('\n window ',i)
    
    maxarea=max(facearea)
    gfaces=[ face for area,face in zip(facearea,facelist) if 
                area>maxarea*ratio]
    return gfaces

def biggestfaces_along_normaldict(wallwindow,wallnormal):
    glassface_bywindowid=defaultdict(list)
    #gpp=GProp_GProps()    
    #print(" wall norm ", wall_norm.Coord())
    for w_id in wallwindow.keys():
        if (w_id in wallnorm.keys()):
            wall_norm=wallnormal[w_id]
        else:
            # window in interior wall
            continue
            
        for win_id in wallwindow[w_id]:
        
            windowshape=create_shape(setting, ifc_file.by_guid(win_id)).geometry
            gfaces=biggestface_along_vector(windowshape,wall_norm)
            glassface_bywindowid[win_id].extend(gfaces)
            
    return glassface_bywindowid
 
def window_in_wall(ifcwall):
    windows=[]
    for op in ifcwall.HasOpenings:
            #print('\n ***',w)
            #print('  ',op)
            for re in op.RelatedOpeningElement.HasFillings:
                #print('Related ', re.RelatedBuildingElement)
                if(re.RelatedBuildingElement.is_a()=='IfcWindow'):
                    windows.append(re.RelatedBuildingElement.id())
    return windows

def link_wall_window(ifcwalls):
    #link window and walls in plain python
    wallwindow=defaultdict(list)
    
    for wall in ifcwalls:
        wallwindow[wall.id()] = window_in_wall(wall)
        
    return wallwindow

class diffuse_mask_on_face:
    def __init__(self,face):
        self._face=face
    
    def compute_mask(self,exposed_building,min_area):
        srf = BRep_Tool().Surface(self._face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
            
              
        t1,p1,areas,mask_vector,vec,sky=compute_diffuse_mask_on_face(exposed_building,self._face,face_norm,min_area)
        #print('\n     diffuse mask ok ')
        
        #self._mask_faces[i].append(shadow_face)
        #self._durations_byfaces[i].append(end-start)
        self._mask_value=mask_vector
        self._masks_areas=areas
        self._directions=vec 
        self._mask_sky=sky
     
    def compute_weighted_mask(self):
        
        m=self._mask_value
        a=self._masks_areas
        ms=self._mask_sky
        
        self._wm=(m*a).sum()/(2*np.pi)
            
        #print(' area sky ', a[ms].sum()/np.pi)
        #print(' area soil ', a[~ms].sum()/np.pi)
            
        self._wm_sky=(m[ms]*a[ms]).sum()/(a[ms].sum())
        self._wm_soil=(m[~ms]*a[~ms]).sum()/(a[~ms].sum())

        print(" wm ", self._wm)
        print(" sky ",self._wm_sky)
        print(" soil ", self._wm_soil)


class diffuse_mask_on_faces:
    """ simple container to hold computation results """
    def __init__(self,lfaces):
        self._lfaces=lfaces
        self._masks=[]
        self._masks_areas=[]
        self._mask_sky=[]
        self._directions=[]
        self._wm=[]
        self._wm_sky=[]
        self._wm_soil=[]                

    def compute_mask(self,exposed_building,min_area):
        for i,gf in enumerate(self._lfaces):
            # re computation of the face normal
            # shoudl be pointing outward
            srf = BRep_Tool().Surface(gf)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(gf.Orientation()==1):
                face_norm.Reverse()
                
                  
            t1,p1,areas,mask_vector,vec,sky=compute_diffuse_mask_on_face(exposed_building,gf,face_norm,min_area)
            print('     diffuse mask ok ')
            
            #self._mask_faces[i].append(shadow_face)
            #self._durations_byfaces[i].append(end-start)
            self._masks.append(mask_vector)
            self._masks_areas.append(areas)
            self._directions.append(vec) 
            self._mask_sky.append(sky)
        
    def compute_weighted_mask(self):
        
        for i,(m,a,d,ms) in enumerate(zip(self._masks,self._masks_areas,self._directions,self._mask_sky)):
            self._wm.append((m*a).sum()/(2*np.pi))
            
            #print(' area sky ', a[ms].sum()/np.pi)
            #print(' area soil ', a[~ms].sum()/np.pi)
            
            self._wm_sky.append((m[ms]*a[ms]).sum()/(a[ms].sum()))
            self._wm_soil.append((m[~ms]*a[~ms]).sum()/(a[~ms].sum()))

        print(" wm ", self._wm)
        print(" sky ",self._wm_sky)
        print(" soil ", self._wm_soil)


class direct_mask_on_face2:
    def __init__(self,face,lsun_dir):
        self._face=face
        self._lsun_dir=lsun_dir
        self._mask_faces=[]
    
    def compute_mask(self,exposed_building,min_area):
        srf = BRep_Tool().Surface(self._face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
            
        for j,sun_dir in enumerate(self._lsun_dir):
            
            mask_face=shadow_caster_ext_on_face(sun_dir,exposed_building,self._face,face_norm,min_area)
            print('\r     sun dir ',j,'/',len(self._lsun_dir)-1,end="",flush=True)
            self._mask_faces.append(mask_face)
    
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



class direct_mask_on_faces:
    """ simple container to hold computation results """
    def __init__(self,lfaces,lsun_dir):
        self._lfaces=lfaces
        self._lsun_dir=lsun_dir
        self._mask_faces=[[] for i in range(len(self._lfaces))]
        self._durations_byfaces=[[]]
        self._mask_area_sum=[]
        self._totalduration=[]
        

    def compute_mask(self,exposed_building,min_area):
        for i,gf in enumerate(self._lfaces):
            # re computation of the face normal
            # shoudl be pointing outward
            srf = BRep_Tool().Surface(gf)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(gf.Orientation()==1):
                face_norm.Reverse()
                
            for j,sun_dir in enumerate(self._lsun_dir):
                start=timer()
                mask_face=shadow_caster_ext_on_face(sun_dir,exposed_building,gf,face_norm,min_area)
                print('\r     sun dir ',j,'/',len(self._lsun_dir)-1,end="",flush=True)
                end=timer()
                self._mask_faces[i].append(mask_face)
                self._durations_byfaces[i].append(end-start)
                
                
        #print(' faces ',self._shadow_faces)
    
    def compute_area_and_ratio(self):
        gpp=GProp_GProps() 
        self._glass_area=0.0
        for gf in self._lfaces :
            brepgprop_SurfaceProperties(gf,gpp)
            self._glass_area+=gpp.Mass()
        
        
        for vector_idx in range(len(self._lsun_dir)):
            area_sum=0.0           
            for face_idx in range(len(self._lfaces)):
                brepgprop_SurfaceProperties(self._mask_faces[face_idx][vector_idx],gpp)
                area_sum+=gpp.Mass()
                
            self._mask_area_sum.append(area_sum)
            #self._totalduration.append( self._durations[face_idx])
        
        # fraction of flux reaching the face
        self._ratio_vector=[1.- a/self._glass_area for a in self._mask_area_sum]
        #print(' shadow area vector ',self._shadow_area_vector)
        #print(' ratio vector ',self._ratio_vector)
        
    def compute_area_and_ratio2(self):
        gpp=GProp_GProps() 
        self._lglass_area=[]
        
        #self._glass_area=0.0
        
        for gf in self._lfaces :
            brepgprop_SurfaceProperties(gf,gpp)
            #self._glass_area+=gpp.Mass()
            self._lglass_area.append(gpp.Mass())
            
        self._lglass_area=np.array(self._lglass_area)
        
        self._shadows_area=np.zeros((len(self._lfaces),len(self._lsun_dir)))
        
        for vector_idx in range(len(self._lsun_dir)):
            for face_idx in range(len(self._lfaces)):
                brepgprop_SurfaceProperties(self._mask_faces[face_idx][vector_idx],gpp)
                self._shadows_area[face_idx,vector_idx]=gpp.Mass()
        self._mask_value = 1.-self._shadows_area/self._lglass_area[:,np.newaxis]        
            
        #self._mask_area_sum_by_face= np.sum(self._shadows_area,axis=1)
        #self._mask_area_sum_by_position = np.sum(self._shadows_area,axis=0)
        
        # fraction of flux reaching the face
        #self._ratio_vector=[1.- a/self._glass_area for a in self._mask_area_sum]
        #print(' shadow area vector ',self._shadow_area_vector)
        #print(' ratio vector ',self._ratio_vector)    
        
    def compute_area_and_ratio_byunion(self):
        """ 
        could be simpler in terms of code but rely on robustness of OCC to compute on more
        complex configurations 
        
        """
        totalface=fuse_listOfShape(self._lfaces)
        gpp=GProp_GProps() 
        brepgprop_SurfaceProperties(totalface,gpp)
        totalarea=gpp.Mass()
        
        ratio=[]
        for vector_idx in range(len(self._lsun_dir)):
            lfaces=[]
            for face_idx in range(len(self._lfaces)):
                f=self._mask_faces[face_idx][vector_idx]
                if not f.IsNull():
                    lfaces.append(f)
            totalshadow=fuse_listOfShape(lfaces)
            brepgprop_SurfaceProperties(totalshadow,gpp)
            shadowarea=gpp.Mass()
            ratio.append(shadowarea/totalarea)
        
        #print(' shadow area vector ',self._shadow_area_vector)
        print(' ratio vector by union',ratio)
        
        
    def compute_complementary_face(self):
        cutter=BOPAlgo_BOP()
        
        self._complementary_faces=[[] for i in range(len(self._lfaces))]
        
        gpp=GProp_GProps() 
        
        #larea=[]
        
        for vector_idx in range(len(self._lsun_dir)):
            #area=0.0           
            for face_idx in range(len(self._lfaces)): 
                mask_face=self._mask_faces[face_idx][vector_idx]
                glass_face=self._lfaces[face_idx]
                
                if not mask_face.IsNull():
                    cutter.Clear()
                    cutter.SetOperation(BOPAlgo_Operation.BOPAlgo_CUT)
                    cutter.AddArgument(glass_face)
                    cutter.AddTool(mask_face)
                    #cutter.SetFuzzyValue(1e-6)
                    cutter.Perform()
                    complementary=cutter.Shape()
                    #print(' cutter ',complementary)
                            
                else :
                    
                    complementary=glass_face
                    
                
                self._complementary_faces[face_idx].append(complementary)
                
               
class shadow_on_faces_byray:
    """ simple container to hold computation results """
    def __init__(self,lfaces,lsun_dir):
        self._lfaces=lfaces
        self._lsun_dir=lsun_dir
        self._shadow_tab=[[] for i in range(len(self._lfaces))]
        self._durations_byfaces=[[]]
        

    def compute_shadow(self,exposed_building,min_area,N):
        
        self._N=N
        for i,gf in enumerate(self._lfaces):
            # re computation of the face normal
            # shoudl be pointing outward
            srf = BRep_Tool().Surface(gf)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(gf.Orientation()==1):
                face_norm.Reverse()
                
            for j,sun_dir in enumerate(self._lsun_dir):
                start=timer()
                #tab,lshape=shadow_caster_ray(sun_dir,exposed_building,gf,face_norm,N)
                tab=shadow_caster_ray(sun_dir,exposed_building,gf,face_norm,N)
                end=timer()
                self._shadow_tab[i].append(tab)
                #print(start,' ',end)
                self._durations_byfaces[i].append(end-start)
                
                
        #print(' faces ',self._shadow_faces)
    
    def compute_area_and_ratio(self):
        
        self._shadow_area_vector=[]
        self._totalduration=[]
        for vector_idx in range(len(self._lsun_dir)):
            area_sum=0.0           
            for face_idx in range(len(self._lfaces)):
                area_sum+=self._shadow_tab[face_idx][vector_idx].sum()
                
            self._shadow_area_vector.append(area_sum)
            #print(self._durations)
        
        
        #self._totalduration.append( self._durations[face_idx])
            
        self._ratio_vector=[ a/(self._N*self._N) for a in self._shadow_area_vector]
        #print(' shadow area vector ',self._shadow_area_vector)
        print(' ratio vector ray',self._ratio_vector)                
        
        



class rtaa_on_faces:
    
    def __init__(self,lfaces,wall_plane,exposed_building,l_sun_dir):
        self._lfaces=lfaces
        self._wall_plane=wall_plane
        #self._original_lfaces=[f for f in lfaces]
        self._exposed_building=exposed_building
        self._lsun_dir=l_sun_dir
        
        
    def cut_exposed_building(self,cut_size,face_ref):
        
        #cut with face_ref as a reference plane of size cut_size
        face=face_ref#self._lfaces[face_ref]
        
        srf = BRep_Tool().Surface(face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(face.Orientation()==1):
            face_norm.Reverse()
                    
        newface= BRepBuilderAPI_MakeFace(plane.Pln(),-cut_size,cut_size,-cut_size,cut_size).Face()
        extrusion = BRepPrimAPI_MakePrism(newface,gp_Vec(face_norm)*10.,False,False).Shape()
        
        intersector=BOPAlgo_BOP()
        intersector.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
        intersector.AddTool(extrusion) 
        intersector.AddArgument(self._exposed_building)
        intersector.Perform()
        return shapes_as_solids([intersector.Shape()])[0]
    
    def adjust_face_to_wall_plane(self,face_ref):
        # translation of the face
        face=face_ref#self._lfaces[face_ref]
        ext_plane=self._wall_plane.Pln()
        srf = BRep_Tool().Surface(face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(face.Orientation()==1):
            face_norm.Reverse()
        
        plane=plane.Pln()
        plane_diff=gp_Vec(ext_plane.Location().XYZ())-gp_Vec(plane.Location().XYZ())
        distance = plane_diff.Dot(gp_Vec(face_norm))
        
        
        translation = gp_Trsf()
        translation.SetTranslation(gp_Vec(face_norm)*distance)
        builder = BRepBuilderAPI_Transform(translation)
        builder.Perform(face)
        #print("Adjsuted face : new face with translation ",distance)
       
        return builder.ModifiedShape(face)
        
    
    def compute_masks(self):
        
        self._ldirect=[]
        self._ldiffuse=[]
        for f in self._lfaces:
         #rtaa_on_faces(exposed_building,sun_to_earth_project)
            cutted_building = self.cut_exposed_building(3.,f)
            adjusted_face = self.adjust_face_to_wall_plane(f)   
            print('direct mask')
            direct_mask = direct_mask_on_face2(adjusted_face,self._lsun_dir)
            direct_mask.compute_mask(cutted_building,1e-3)
            direct_mask.shadow_areas()
            direct_mask.face_area()
            direct_mask.mask_ratio()
            print('\n')
            
            print('diffuse mask ')
            diffuse_mask=diffuse_mask_on_face(adjusted_face)
            diffuse_mask.compute_mask(cutted_building,1e-3)
            diffuse_mask.compute_weighted_mask()
            print('\n')
            
            self._ldirect.append(direct_mask)
            self._ldiffuse.append(diffuse_mask)
            
        
    def compute_cm(self,irradiance):
        
        albedo=.2
        origin = gp_Pnt(0.0,0.0,0.0)
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
            
            directflux = [face_norm.Dot(s.Reversed()) for s in self._lsun_dir]
                        
            Drp = irradiance['DNI']*[f if f>-1.e-5 else 0.0 for f in directflux]
            Dfp = irradiance['DHI']*(1+cosbeta)*.5
            Rrp = irradiance['GHI']*albedo*(1-cosbeta)*.5
            
            mDrp = Drp*dir._mask_ratio
            mDfp = Dfp*diff._wm_sky
            mRrp = Rrp*diff._wm_soil
            
            face_irr=dict()
            face_irr['Drp']=Drp
            face_irr['Dfp']=Dfp
            face_irr['Rrp']=Rrp
            
            face_irr['mDrp']=mDrp
            face_irr['mDfp']=mDfp
            face_irr['mRrp']=mRrp
            
            face_irr['mask']=dir._mask_ratio
            
            self._ldf_irr.append(pd.DataFrame(face_irr))
            
            masked_irr   = mDrp + mDfp + mRrp 
            unmasked_irr =  Drp +  Dfp + Rrp
            
            self._cm_byface.append( masked_irr.sum()/unmasked_irr.sum() )
            
        
        all_values=pd.concat(self._ldf_irr,axis=0)
        total_irr =all_values[['Drp','Dfp','Rrp']].sum().sum()
        total_m_irr =all_values[['mDrp','mDfp','mRrp']].sum().sum()
        
        self._cm= total_m_irr/total_irr
        
    
        
        
class project_location:
    def __init__(self):
        pass
        
    def set_northing_from_ifc(self,ifc_file):
        repr_context = ifc_file.by_type('IfcGeometricRepresentationContext',False)
        project_repre = repr_context[0]
        true_north = project_repre.TrueNorth
        tn_X,tn_Y= true_north.DirectionRatios
        
        # true north vector in project coordinate system
        tn_vec = gp_Vec(tn_X,tn_Y,0.0)
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Xaxis = gp_Ax1(origin,gp_Dir(1.0,0.0,0.0))
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        # transformation to apply to convert in project coordinates
        # any vector expressed in world coordinate (sun direction)
        #self._tn_angle_sgn = self._tn_vec.AngleWithRef(gp_Vec(Yaxis.Direction()),gp_Vec(Zaxis.Direction()))
        self._tn_angle     = tn_vec.Angle(gp_Vec(Yaxis.Direction()))
        
        print("Angle true North : ",tn_angle)
        
        #print("Signed angle : ", self._tn_angle_sgn)
        """
        self._tn_rotation = gp_Trsf()
        self._tn_rotation.SetRotation(Zaxis,self._tn_angle)
        print("TN vector : ",(Yaxis.Direction().Transformed(self._tn_rotation).Coord()))
        """
        
    def update_northing_from_angle(self,new_tn):
        #self._tn_angle_sgn = new_tn
        self._tn_angle     = np.abs(new_tn)
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        #print("
        """
        self._tn_rotation = gp_Trsf()
        self._tn_rotation.SetRotation(Zaxis,self._tn_angle)
        self._tn_vec = gp_Vec(Yaxis.Direction()).Transformed(self._tn_rotation)
        print("Updated angle true North : ",self._tn_angle)
        #print("Updated Signed angle : ", self._tn_angle_sgn)
        """
        
    
    def set_location_from_ifc(self,ifc_file):
        ## Site location
        ifcsite = ifc_file.by_type('IfcSite')[0]
        h,m,s,ms = ifcsite.RefLatitude
        self._latitude = h+m/60+(s+ms*1e-6)/3600
        h,m,s,ms = ifcsite.RefLongitude
        self._longitude= h+m/60+(s+ms*1e-6)/3600
        print("latitude : ",self._latitude)
        print("longitude: ",self._longitude)
        
        ## datetime to compute shadow
        tf=tf = TimezoneFinder()
        self._tz = tf.timezone_at(lng=self._longitude, lat=self._latitude) 
        print("TimeZone : ",self._tz)

        
        
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
        
    
    
    

if __name__ == "__main__":

    # due to some bugs in ipython parsing
    __import__("logging").getLogger("parso.python.diff").setLevel("INFO")
    __import__("logging").getLogger("parso.cache").setLevel("INFO")
    __import__("logging").getLogger("asyncio").setLevel("INFO")

    

    # Initialize a graphical display window (from ifcos)

    setting=ifcopenshell.geom.settings()
    setting.set(setting.USE_PYTHON_OPENCASCADE, True)

    filename = 'data/Rtaa_validation_run_onewindow.ifc'
    filename = 'data/Rtaa_validation_run.ifc'
    #ifc_file= ifcopenshell.open('data/simple_reunion_northaligned.ifc')
    #ifc_file= ifcopenshell.open('data/Rtaa_validation_run.ifc')

    ifc_file= ifcopenshell.open(filename)
    
    tag2ratio={}
    if filename == 'data/Rtaa_validation_run.ifc':
        tag2ratio['257017']=0.1
        tag2ratio['273402']=0.2
        tag2ratio['273433']=0.3
        tag2ratio['273468']=0.4
        tag2ratio['273497']=0.5
        tag2ratio['273528']=0.6
        tag2ratio['273603']=0.7
        tag2ratio['273636']=0.8
        tag2ratio['273680']=0.9
        tag2ratio['273718']=1.0
        

    ifcwalls=ifc_file.by_type('IfcWall')
    ifcspaces=ifc_file.by_type('IfcSpace')
    ifcwindows=ifc_file.by_type('IfcWindow')
    ifcslabs=ifc_file.by_type('IfcSlab')
    ifcproxys=ifc_file.by_type('IfcBuildingElementProxy')
    
        
    # partial building to compute external shell and exterior wall
    wall_shapes  = [create_shape(setting, x).geometry for x in ifcwalls if x.Representation is not None]
    space_shapes = [create_shape(setting, x).geometry for x in ifcspaces if x.Representation is not None]
    core_shapes  = wall_shapes+space_shapes
    core_solids  = shapes_as_solids(core_shapes)
        
    # complete building to compute shadow on
    ifcextension= []+ifcslabs+ifcproxys
    extension_shapes = [create_shape(setting, x).geometry for x in ifcextension if x.Representation is not None]
    extension_solids =  shapes_as_solids(extension_shapes)

    building_shapes= core_shapes + extension_shapes
    building_solids= core_solids + extension_solids
    exposed_building = fuse_listOfShape(building_solids)
    exposed_building = shapes_as_solids([exposed_building])[0]
    
    external_shell= get_external_shell(core_solids)
    
    windows_by_wall_id = dict()
    for wall in ifcwalls:
        windows_by_wall_id[wall.id()] = window_in_wall(wall)

    # will only contain wall id that considered as exterior wall of the building
    # if id is not in the keys, the wall could be considered as interior
    normal_by_exterior_wall_id = dict()
    plane_by_exterior_wall_id=dict()
    for (w_id,ws) in zip(windows_by_wall_id.keys(),wall_shapes):
        wall_norm,plane_ext=exterior_wall_normal_and_plane(ws,external_shell)
        normal_by_exterior_wall_id[w_id]=wall_norm
        plane_by_exterior_wall_id[w_id]=plane_ext
    
    
    glassface_bywindowid=defaultdict(list)
    glassface_extplane=dict()
    for w_id in windows_by_wall_id.keys():
        if (w_id in normal_by_exterior_wall_id.keys()):
            wall_norm=normal_by_exterior_wall_id[w_id]
        else:
            # window in interior wall
            continue
            
        for win_id in windows_by_wall_id[w_id]:
        
            windowshape=create_shape(setting, ifc_file.by_guid(win_id)).geometry
            gfaces=biggestface_along_vector(windowshape,wall_norm)
            glassface_bywindowid[win_id].extend(gfaces)
            glassface_extplane[win_id]=plane_by_exterior_wall_id[w_id]
        
    
    irradiance=pd.read_csv('data/irradiance_data/2659942_-21.15_55.62_2018.csv',skiprows=[0,1])
    # creating datetime object
    col= irradiance.columns
    dt=pd.to_datetime(irradiance[col[:5]])
    # parsing usefull datas
    irr_col=['Clearsky DHI','Clearsky DNI', 'Clearsky GHI',
            'DHI', 'DNI','GHI']
    irradiance_only=irradiance[irr_col]
    irradiance_only=irradiance_only.assign(time=dt.values)
    
    north_mask = (dt<'2018-04-30')*(dt>'2018-02-01')
    eso_mask = (dt<'2018-02-28')*(dt>'2018-01-01')
    
        
    tn_angle=[0.0,np.pi*.5,np.pi,3.*np.pi*.5]
    
    orientations_name=['NORD','EST','SUD','OUEST']
    time_mask = [north_mask,eso_mask,eso_mask,eso_mask]
    
    tn_angle=[0.0]
    proj_loc=project_location()
    proj_loc.set_location_from_ifc(ifc_file)
    proj_loc.set_northing_from_ifc(ifc_file)
    
    res=[]
    for tn,mask in zip(tn_angle,time_mask):
        lcm=[]
        for (k,win_id) in enumerate(glassface_bywindowid.keys()):
            
            lglassfaces=glassface_bywindowid[win_id]
            print(' window id ', win_id)
            #if win_id!=626:
            #    continue
            
            
            sample= irradiance_only[mask]
            sample= sample[:100]
            
            dr = pd.DatetimeIndex(sample['time']) 
            
            #proj_loc.update_northing_from_angle(tn)
            sun_to_earth_project,zen_vec,az_vec=proj_loc.sun_vectors(dr)
                
            rtaa=rtaa_on_faces(lglassfaces,glassface_extplane[win_id],exposed_building,sun_to_earth_project)
            rtaa.compute_masks()
            rtaa.compute_cm(sample)
            
            lcm.append(rtaa._cm)
        res.append(lcm)
    
    
    rtaa_data= pd.read_excel('data/h85_l50_ec15.xlsx')
    
    x=np.arange(0.1,1.1,.1)
    colors=['r','g','b','k']
    for data,name,c in zip(res,orientations_name,colors):
        plt.plot(x,data,'-',color=c,lw=2,label=name+'_computed')
        plt.plot(x,rtaa_data[name],'--',color=c,lw=2,label=name+'_ref')
    
    plt.legend()
    plt.show()
 
        
    
    """
    fig = plt.figure()#figsize=(4., 4.))
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(1,len(ldiff)),  # creates 2x2 grid of axes
                 axes_pad=0.1,  # pad between axes in inch.
                 )
    #lmask=[ x[-1] for x in ldiff]
    skyonly=False
    for ax,(dir,sky,mask) in zip(grid, ldiff):
        
        skymask=sky[0]
        x=dir[0][:,0]
        y=dir[0][:,1]
        z=(mask[0]).flatten()
        
        if skyonly:
            x=x[skymask]
            y=y[skymask]
            z=z[skymask]
        
        triang = tri.Triangulation(x, y)
        ax.scatter(x,y,s=z)
        ax.set_aspect('equal')
        tcf = ax.tricontourf(triang, z,cmap=plt.colormaps.get_cmap('viridis'))
        ax.scatter(x,y)

    fig.colorbar(tcf)
    plt.show()
    
    # valid for the file test
    x=np.arange(0.1,1.1,.1)
    rta_north= [0.95,0.88,0.83,0.8,0.77,0.76,.75,.75,0.74,0.74]
    plt.plot(x,rta_north,'-o',label='rta')
    plt.plot(x,lcm,'-o',label='ifc')
    plt.legend()
    plt.show()
    """
    
    """
    result=pd.DataFrame()
    result.index = df.index
    for id,sof in sofdict.items():
        name=ratio_by_id[id]
        result[id]=sof._ratio_vector
    
    plt.plot(result,'o-')
    plt.show()
    """
    """
    def rgb_color(r, g, b):
        return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
    x=50/256
    gray=rgb_color(x, x, x)

    display, start_display, add_menu, add_function_to_menu = init_display()

    display.DisplayShape(exposed_building,color=gray,transparency=0.9)
    #[display.DisplayShape(s,color=gray,transparency=0.9) for s in building_shapes]
    
    #[display.DisplayShape(s,color='BLUE',transparency=0.1) for s in lglass]
    
    #lvec=[gp_Vec(v) for v in earth_to_sun_project]
    
    [display.DisplayVector(v*3.,origin) for v in lvector_diff]
    #display.DisplayVector(tn_vec*2.,origin)
    #display.DisplayVector(tn_proj*4.,origin)
    #display.DisplayShape(external_shell,color='RED',transparency=0.5)
    
    #for sof in lsof:
    #   [display.DisplayShape(s,transparency=0.1,color='BLACK') for s in sof._shadow_faces]
    #   [display.DisplayShape(s,transparency=0.1,color='YELLOW') for s in sof._complementary_faces]
       
    #[display.DisplayShape(s,color='RED',transparency=0.9) for s in lext2[:]] 
    #[display.DisplayShape(s,color='BLUE',transparency=0.5) for s in linter[::3]] 

    #[display.DisplayShape(s,color='GREEN',transparency=0.1) for s in lshells] 

    #[display.DisplayShape(s,color='BLUE',transparency=0.1) for s in lshapee if s is not None] 
    #[display.DisplayShape(s,color='RED',transparency=0.0) for s in lcyl]

    #s=lsof[0]._complementary_faces[0]

    display.FitAll()
    #ifcopenshell.geom.utils.main_loop()
    start_display()
    """
    
    

