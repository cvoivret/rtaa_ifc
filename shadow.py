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
import gc

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
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Copy
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
lmod=[]
lntry=[]
lprec=[]
def shadow_caster_ext(sun_dir,mybuilding,myface,theface_norm,min_area = 1e-3):
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

    
    if sun_dir.Z()>0.:
        #print('Z positive')
        return theface
    
    
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
        
    lntry.append(ntry)    
    lprec.append(precision)
        
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

def exterior_wall_normal(wall_shape,external_shell):
    
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
        
        return wall_norm
    

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


class shadow_on_faces:
    """ simple container to hold computation results """
    def __init__(self,lfaces,lsun_dir):
        self._lfaces=lfaces
        self._lsun_dir=lsun_dir
        self._shadow_faces=[[] for i in range(len(self._lfaces))]
        self._durations_byfaces=[[]]
        

    def compute_shadow(self,exposed_building,min_area):
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
                shadow_face=shadow_caster_ext(sun_dir,exposed_building,gf,face_norm,1.e-3)
                print('     sun dir ',j,'/',len(self._lsun_dir))
                end=timer()
                self._shadow_faces[i].append(shadow_face)
                self._durations_byfaces[i].append(end-start)
                
                
        #print(' faces ',self._shadow_faces)
    
    def compute_area_and_ratio(self):
        gpp=GProp_GProps() 
        self._glass_area=0.0
        for gf in self._lfaces :
            brepgprop_SurfaceProperties(gf,gpp)
            self._glass_area+=gpp.Mass()
        
        self._shadow_area_vector=[]
        self._totalduration=[]
        for vector_idx in range(len(self._lsun_dir)):
            area_sum=0.0           
            for face_idx in range(len(self._lfaces)):
                brepgprop_SurfaceProperties(self._shadow_faces[face_idx][vector_idx],gpp)
                area_sum+=gpp.Mass()
                
            self._shadow_area_vector.append(area_sum)
            #self._totalduration.append( self._durations[face_idx])
            
        self._ratio_vector=[ a/self._glass_area for a in self._shadow_area_vector]
        #print(' shadow area vector ',self._shadow_area_vector)
        print(' ratio vector ',self._ratio_vector)
        
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
                f=self._shadow_faces[face_idx][vector_idx]
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
                shadow_face=self._shadow_faces[face_idx][vector_idx]
                glass_face=self._lfaces[face_idx]
                
                if not shadow_face.IsNull():
                    cutter.Clear()
                    cutter.SetOperation(BOPAlgo_Operation.BOPAlgo_CUT)
                    cutter.AddArgument(glass_face)
                    cutter.AddTool(shadow_face)
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
    
    ratio_by_id={w.id():tag2ratio[w.Tag]for w in ifcwindows}
    
    """
    doors=ifc_file.by_type('IfcDoor')
    opening=ifc_file.by_type('IfcOpeningElement')
    storeys=ifc_file.by_type('IfcBuildingStorey')
    roof=ifc_file.by_type('IfcRoof')
    """
    
    repr_context = ifc_file.by_type('IfcGeometricRepresentationContext',False)
    project_repre = repr_context[0]
    true_north = project_repre.TrueNorth
    tn_X,tn_Y= true_north.DirectionRatios
    tn_vec = gp_Vec(tn_X,tn_Y,0.0)
    
    origin = gp_Pnt(0.0,0.0,0.0)
    Xaxis = gp_Ax1(origin,gp_Dir(1.0,0.0,0.0))
    Yaxis = gp_Ax1(origin,gp_Dir(0.0,1.0,0.0))
    Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
    
    
    # transformation to apply to convert in project coordinates
    # any vector expressed in world coordinate (sun direction)
    angle = tn_vec.Angle(gp_Vec(Yaxis.Direction()))
    world_to_project = gp_Trsf()
    world_to_project.SetRotation(Zaxis,-angle)
    print("Angle true North : ",angle)
    tn_proj = tn_vec.Transformed(world_to_project)
    
    ## Site location
    ifcsite = ifc_file.by_type('IfcSite')[0]
    h,m,s,ms = ifcsite.RefLatitude
    latitude = h+m/60+(s+ms*1e-6)/3600
    h,m,s,ms = ifcsite.RefLongitude
    longitude= h+m/60+(s+ms*1e-6)/3600
    print("latitude : ",latitude)
    print("longitude: ",longitude)
    
    ## datetime to compute shadow
    tf=tf = TimezoneFinder()
    tz = tf.timezone_at(lng=longitude, lat=latitude)  # 'Europe/Berlin'    
    print("TimeZone : ",tz)
    
    """
    days=pd.date_range(start='15/1/2020',end='15/12/2020',freq='M')
    hoursofdays =[ pd.date_range(start=d,periods=24,freq='H') for d in days]
    dr=hoursofdays[0]
    for hd in hoursofdays[1:]:
        dr=dr.append( hd)
    """
    #dr = pd.date_range(start='2020/5/25',end='2020/6/7',freq='H',inclusive="neither")
    dr = pd.date_range(start='2020/1/1',end='2020/12/31',freq='H',inclusive="neither")
    #dr=dr[1880:2000]

    dr_proj = dr.tz_localize(tz)
    dr_proj_utc = dr_proj.tz_convert("UTC")
    
    az_vec,zen_vec=sunpos.sunpos(dr_proj_utc,latitude,longitude,0)[:2]
    elev_vec=90-zen_vec     
    
    #-21.34053399695468, 55.49058057798694
    
    
    # take X axis

    
    lvector=[]
    earth_to_sun_project=[]
    sun_to_earth_project=[]
    # create transform along Xaxis based on altitude
    for zen,az in zip(zen_vec,az_vec):
        RotX = gp_Trsf()
        RotX.SetRotation(Xaxis,np.deg2rad(90-zen))
        elev_dir = Yaxis.Transformed(RotX)
        # create rotation around Z axis for azimuth
        RotZ=gp_Trsf()
        RotZ.SetRotation(Zaxis,np.deg2rad(-az))
        sun_axis=elev_dir.Transformed(RotZ)
        # need to take into account True north
        sun_direction=sun_axis.Direction()
        lvector.append(sun_direction.Coord())
        earth_to_sun_project.append(sun_direction)
        sun_to_earth_project.append(sun_direction.Reversed())
    
    sun_vec_world=np.array(lvector)
    df=pd.DataFrame(data=np.column_stack([az_vec,elev_vec,sun_vec_world]),
                    index=dr_proj,
                    columns=['azimuth','elevation','Vx','Vy','Vz'])
    
    
    #sun_to_earth_project = [ d.Reversed() for d in earth_to_sun_project]
    
    
    
    #df['local_datetime']=dr_proj
    
    #df_daylight = df[ df.elevation>0.0]
    
    #df_daylight=df_daylight[::1]
    
    #earth_to_sun_project = [gp_Dir(*vec).Transformed(world_to_project)
     #                   for vec in df.values[:,2:5]]
    #earth_to_sun_project2 = [gp_Vec(vec) for vec in earth_to_sun_project]
    
    #sun_to_earth_project = [ d.Reversed() for d in earth_to_sun_project]
    

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
    for (w_id,ws) in zip(windows_by_wall_id.keys(),wall_shapes):
        normal_by_exterior_wall_id[w_id]=exterior_wall_normal(ws,external_shell)
    
    
    glassface_bywindowid=defaultdict(list)
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
    #lsof=[]
    sofdict=dict()
    #lhalf=[]
    #lglass=[]
    #lext=[]
    linter=[]
    for (k,win_id) in enumerate(glassface_bywindowid.keys()):
        lglassfaces=glassface_bywindowid[win_id]
        print(' window id ', win_id)
        #if win_id!=1112:
        #    continue
        face= lglassfaces[0]
        #lglass.append(face)
        
        srf = BRep_Tool().Surface(face)
        
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(face.Orientation()==1):
            face_norm.Reverse()
        
        umin,umax,vmin,vmax=breptools_UVBounds(face)
        centerXYZ = srf.Value(0.5*(umax-umin),0.5*(vmax-vmin))
        #print(centerXYZ.Coord())
        centerXYZ.Translate(gp_Vec(face_norm)*10)
        size=3.
        newface= BRepBuilderAPI_MakeFace(plane.Pln(),-size,size,-size,size).Face()
        halfspace = BRepPrimAPI_MakeHalfSpace(newface,centerXYZ).Solid()
        #lhalf.append(halfspace)
        extrusion = BRepPrimAPI_MakePrism(newface,gp_Vec(face_norm)*10.,False,False).Shape()
        #lext.append(extrusion)
        
        intersector=BOPAlgo_BOP()
        intersector.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
        intersector.AddTool(extrusion) 
        intersector.AddArgument(exposed_building)
        intersector.Perform()
        intersection=shapes_as_solids([intersector.Shape()])[0]
        
        print(" inter mod init",intersection.Modified())
        """
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)

        display, start_display, add_menu, add_function_to_menu = init_display()
        #[display.DisplayShape(s,color=gray,transparency=0.5) for s in building_shapes]
        [display.DisplayShape(s,color='BLUE',transparency=0.0) for s in lglassfaces]


        display.DisplayShape(intersection,color='RED',transparency=0.5)
        display.FitAll()
        start_display()
        """
        
        sof=shadow_on_faces(lglassfaces,sun_to_earth_project)
        sof.compute_shadow(intersection,1e-3)
        sof.compute_area_and_ratio()
        sof.compute_complementary_face()
        sofdict[win_id]=sof
        #lsof.append(sof)
    
    
    result=pd.DataFrame()
    result.index = df.index
    for id,sof in sofdict.items():
        name=ratio_by_id[id]
        result[id]=sof._ratio_vector
    
    plt.plot(result,'o-')
    plt.show()
    """
    def rgb_color(r, g, b):
        return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
    x=50/256
    gray=rgb_color(x, x, x)

    display, start_display, add_menu, add_function_to_menu = init_display()

    display.DisplayShape(exposed_building,color=gray,transparency=0.9)
    #[display.DisplayShape(s,color=gray,transparency=0.9) for s in building_shapes]
    
    [display.DisplayShape(s,color='BLUE',transparency=0.1) for s in lglass]
    
    lvec=[gp_Vec(v) for v in earth_to_sun_project]
    
    [display.DisplayVector(v*2.,origin) for v in lvec]
    #display.DisplayVector(tn_vec*2.,origin)
    #display.DisplayVector(tn_proj*4.,origin)
    #display.DisplayShape(external_shell,color='RED',transparency=0.5)
    
    #for sof in lsof:
    #   [display.DisplayShape(s,transparency=0.1,color='BLACK') for s in sof._shadow_faces]
    #   [display.DisplayShape(s,transparency=0.1,color='YELLOW') for s in sof._complementary_faces]
       
    [display.DisplayShape(s,color='RED',transparency=0.9) for s in lext2[:]] 
    #[display.DisplayShape(s,color='BLUE',transparency=0.5) for s in linter[::3]] 

    #[display.DisplayShape(s,color='GREEN',transparency=0.1) for s in lshells] 

    #[display.DisplayShape(s,color='BLUE',transparency=0.1) for s in lshapee if s is not None] 
    #[display.DisplayShape(s,color='RED',transparency=0.0) for s in lcyl]

    #s=lsof[0]._complementary_faces[0]

    display.FitAll()
    #ifcopenshell.geom.utils.main_loop()
    start_display()
    """
    
    

