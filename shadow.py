from collections import defaultdict
import itertools
import os
from itertools import compress,tee,combinations
from collections import Counter
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
from OCC.Core.gp import gp_Pnt2d,gp_Dir2d,gp_Vec2d,gp_Lin2d

from OCC.Core.Geom import Geom_Plane,Geom_Line
from OCC.Core.GeomAPI import GeomAPI_IntCS,GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomProjLib import geomprojlib_Project
from OCC.Core.GeomProjLib import geomprojlib_ProjectOnPlane


from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.TopTools import TopTools_ListOfShape,TopTools_IndexedMapOfShape
from OCC.Core.TopExp import topexp_MapShapes
from OCC.Core.TopAbs import TopAbs_SOLID,TopAbs_FACE,TopAbs_SHELL,TopAbs_WIRE

from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer

from OCC.Core.GeomAdaptor import GeomAdaptor_Curve


from OCC.Core.ProjLib import ProjLib_ProjectOnPlane
from OCC.Core.ProjLib import projlib_Project

 
from OCC.Core.BOPAlgo import BOPAlgo_BOP,BOPAlgo_Operation
from OCC.Core.BOPAlgo import BOPAlgo_CellsBuilder
from OCC.Core.BOPTools import BOPTools_AlgoTools_OrientFacesOnShell

from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib

import OCC.Core.ShapeFix as ShapeFix_Shape
from OCC.Core.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve
from OCC.Core.GeomAPI import geomapi_To3d
from OCC.Core.Geom2d import Geom2d_Line
from OCC.Core.ElCLib import elclib_Parameter,elclib_To3d,elclib_Value

from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector

from OCC.Core.Standard import standard_Purge

from OCC.Extend.DataExchange import write_stl_file

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


setting=ifcopenshell.geom.settings()
setting.set(setting.USE_PYTHON_OPENCASCADE, True)

    
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


    
def ifcelement_as_solid(element):
    if element.Representation is None:
        return None
    shape=create_shape(setting, element).geometry
    #print(shape)
    solids=shapes_as_solids([shape])
    #print(solids)
    solid=None
    if( len(solids)==1):
        solid = solids[0]
    else :
        solid = fuse_listOfShape(solids)
    
    shu=ShapeUpgrade_UnifySameDomain(solid)
    shu.Build()
    solid=shu.Shape()    
    #print(solid)
    
    # add unifiysamedomain ?
    
    return solid
        
        
     


def get_external_shell(lshape):
    """
    try to identigy a shell (set of face) that limit the inside and the outside of the building.
    Basically, wall and room must be part of the solid list input.
    Build the boundign box of the whole model and enlarge it.
    
    """
       
    #unionize solids
    unionsolid=fuse_listOfShape(lshape)
    shu=ShapeUpgrade_UnifySameDomain(unionsolid)
    shu.Build()
    unionsolid=shu.Shape()
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
    
    shu=ShapeUpgrade_UnifySameDomain(commonshell)
    shu.Build()
    commonshell=shu.Shape()
    
    return commonshell


def compute_direct_mask_on_face(sun_dir,mybuilding,myface,theface_norm,min_area = 1e-3):
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
    #print(' DOT ' , theface_norm.Dot(sun_dir),flush=True)
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
            
        
        
        #print(" precision update ", precision,flush=True)
        
    #lntry.append(ntry)    
    #lprec.append(precision)
        
    intersection=cb.Shape()
    
    intersection_faces=list(TopologyExplorer(intersection).faces())
    #print('N iter ',len(intersection_faces),flush=True)
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
            #if(ff.Orientation()==1):
            #    ff.Reverse()
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
        lmask_sky.append( final.Dot(gp_Vec(Zaxis.Direction()))>=0.)
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
    
    
    # compute ratio for each direction
    dmof=direct_mask_on_face(theface,lvec)
    dmof.compute_mask(mybuilding,min_area)
    dmof.shadow_areas()
    dmof.face_area()
    dmof.mask_ratio()
    
    mask_vector = np.array(dmof._mask_ratio)
    mask_sky = np.array(lmask_sky)
    
    return t1,p1,areas,mask_vector,projected_vec,mask_sky,lvec 
    
    # could use sof with discretization vector instead of sun_dir
    # Sum up and return




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
    
    #print(commonshell2)
    
    faces=list(TopologyExplorer(commonshell2).faces())
    #print(len(faces))
    # exteriro wall !!
    if len(faces)>0:
        
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
        #print(norm_map)
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
    #print(" nb face par fenetre ", len(faces))
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
    if(len(facearea)==0):
        # no faces along given vector
        return []
    
    totalarea=sum(facearea)
    
    count=Counter(facearea)
    #area_by_count={area: c*area for area,c in count.items()}
    area_by_count2={area: c*area/totalarea for area,c in count.items()}
    
    significative_area = [area for area,v in area_by_count2.items() if v >.2]
    gfaces = [f for a,f in zip(facearea,facelist) if a in significative_area]
    """
    print(count)
    print(area_by_count)
    print(area_by_count2)
    print(significative_area)
    print([a for a in facearea if a in significative_area])
    """
    # une grande face 
    
     
    
    #maxarea=max(facearea)
    #facearea.sort()
    #print(" face area sorted ",facearea)
    #gfaces=[ face for area,face in zip(facearea,facelist) if 
    #           area>maxarea*ratio]
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
 
def opening_in_wall(ifcwall,take_window=True,take_doors=True):
    opening_filler=[]
    #print('\n ***',ifcwall)
    for op in ifcwall.HasOpenings:
        #print('   opening : ',op)
        for re in op.RelatedOpeningElement.HasFillings:
            #print(' Filling ', re.RelatedBuildingElement.is_a())
            if(take_window & (re.RelatedBuildingElement.is_a()=='IfcWindow')):
                opening_filler.append(re.RelatedBuildingElement.id())
            elif (take_doors & (re.RelatedBuildingElement.is_a()=='IfcDoor')):
                opening_filler.append(re.RelatedBuildingElement.id())
    return opening_filler


class diffuse_mask_on_face:
    def __init__(self,face):
        self._face=face
    
    def compute_mask(self,exposed_building,min_area):
        srf = BRep_Tool().Surface(self._face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(self._face.Orientation()==1):
            face_norm.Reverse()
            
              
        t1,p1,areas,mask_vector,vec,sky,lvec=compute_diffuse_mask_on_face(exposed_building,self._face,face_norm,min_area)
        #print('\n     diffuse mask ok ')
        
        #self._mask_faces[i].append(shadow_face)
        #self._durations_byfaces[i].append(end-start)
        self._mask_value=mask_vector
        self._masks_areas=areas
        self._directions=vec 
        self._mask_sky=sky
        self._lvec = lvec
        self._ldot = -1.*np.array([ face_norm.Dot(v) for v in lvec])
        
        
        #self.display(exposed_building)
    

    def compute_weighted_mask(self):
        
        m=self._mask_value
        a=self._masks_areas
        ms=self._mask_sky
        dot=self._ldot
        
        self._wm=(m*dot*a).sum()/(a*dot).sum()
        self._wm_sky=(m[ms]*a[ms]*dot[ms]).sum()/((a[ms]*dot[ms]).sum())
        self._wm_soil=(m[~ms]*a[~ms]*dot[~ms]).sum()/((a[~ms]*dot[~ms]).sum())

        #print(m[ms])
        #print(a[ms])
        print("\n wm ", self._wm)
        print(" sky ",self._wm_sky)
        print(" soil ", self._wm_soil)
        
        
    def display(self,exposed_building):
        
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
        
        lvec=[gp_Vec(d) for d in self._lvec]
        # display of sun vectors around the face
        pos =[ center.Translated(v*(-3.)) for v in lvec]
        
        
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.DisplayVector(face_norm,center)
        [display.DisplayVector(v,p) for v,p in zip(lvec,pos)]
        display.DisplayShape(exposed_building,color=gray,transparency=0.5)
        display.DisplayShape(self._face,color='BLUE',transparency=0.5)
        
        display.FitAll()
        start_display()

      
        


class direct_mask_on_face:
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
            
        #self.display(exposed_building)
        #csdc
        
        for j,sun_dir in enumerate(self._lsun_dir):
            # sun under the horizon line (nighttime)
            #if sun_dir.Coord()[2]>0.0 :
            #    mask_face=self._face
            #else :
            mask_face=compute_direct_mask_on_face(sun_dir,exposed_building,self._face,face_norm,min_area)
            
            print('\r     sun dir ',j,'/',len(self._lsun_dir)-1,end="",flush=True)
            self._mask_faces.append(mask_face)
        
    def display(self,exposed_building):
        
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


# the list of faces correspond to a window
class rtaa_on_faces:
    
    def __init__(self,lfaces,wall_plane,exposed_building,l_sun_dir,cut=True,adjust=True):
        self._lfaces=lfaces
        self._wall_plane=wall_plane
        #self._original_lfaces=[f for f in lfaces]
        self._exposed_building=exposed_building
        self._lsun_dir=l_sun_dir
        self._cut = cut
        self._adjust = adjust
        
    
    
    
    def cut_exposed_building(self,cut_size,face_ref):
        
        #cut with face_ref as a reference plane of size cut_size
        face=face_ref#self._lfaces[face_ref]
        
        srf = BRep_Tool().Surface(face)
        plane = Geom_Plane.DownCast(srf)
        face_norm = plane.Axis().Direction()
        if(face.Orientation()==1):
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
            direct_mask.compute_mask(cutted_building,1e-3)
            direct_mask.shadow_areas()
            direct_mask.face_area()
            direct_mask.mask_ratio()
            #direct_mask.display(cutted_building)
            print('\n')
            
            print('diffuse mask ')
            diffuse_mask=diffuse_mask_on_face(adjusted_face)
            diffuse_mask.compute_mask(cutted_building,1e-3)
            diffuse_mask.compute_weighted_mask()
            print('\n')
            
            self._ldirect.append(direct_mask)
            self._ldiffuse.append(diffuse_mask)
            
            self._lcuted.append(cutted_building)
            self._ladjusted.append(adjusted_face)
    
    def display(self):
        
        
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
        # direct irradiance
        
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
            #print(np.arccos(directflux))
            #print(Rb)
            #print(irradiance['DNI'])
            #maskratio[not_exposed]=0.0
            
            #actual_ghi = irradiance['DNI']*np.sin(np.pi-zenith)+irradiance['DHI']
            #actual_ghi[not_exposed]=irradiance['DHI'][not_exposed]
            
            Drp = data_irradiance['DNI']*directflux
            Dfp = data_irradiance['DHI']*(1+cosbeta)*.5
            Rrp = data_irradiance['GHI']*albedo*(1-cosbeta)*.5
            # RTAA rule
            Drp[not_exposed]=0.0
            #Dfp[not_exposed]=0.0
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
        relative_path = "data/meteo_rtaa.xlsx"
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
    
   
    
    
class rtaa_solar_study:
    
    def __init__(self,ifcfilename):
        setting=ifcopenshell.geom.settings()
        setting.set(setting.USE_PYTHON_OPENCASCADE, True)
        self._ifc_file= ifcopenshell.open(ifcfilename)
        self._solar_elements=dict()
        self._building_elements=dict()
        print(" ifc file ",self._ifc_file)
        
        self._proj_loc=project_location()
        self._proj_loc.set_location_from_ifc(self._ifc_file)
        self._proj_loc.set_northing_from_ifc(self._ifc_file)
        self._proj_loc.load_irradiance()
    
    def config_from_file(json_file):
        # read ids and class to set elements of interest and buildings
        # add the element 
        #build the geometry
        pass
    
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
                print("Error, solar element not wall window or door")
            self._hosting_solar_id[id]=hosting_wall_id
            
            #print(" wall norm plane ",hosting_wall_id," ",self._wall_norm_plane[hosting_wall_id])
            if(self._wall_norm_plane[hosting_wall_id]):
                wall_norm,plane = self._wall_norm_plane[hosting_wall_id]
                #print(wall_norm.Coord())
                #print(id)
            
                faces = biggestface_along_vector(shape,wall_norm)
                if len(faces)>0:
                    self._solar_faces[id]=faces
            else :
                print(" No normal found for this window ",id)
        #print(self._solar_faces)  

        
    def _building_solid(self):
        
        building_shapes = { x.id():ifcelement_as_solid(x) 
                            for x in self._building_elements.values() }
        
        self._building= fuse_listOfShape(list(building_shapes.values()))
        # unify ?
        
        # get spaces geometries
        ifcspaces=self._ifc_file.by_type('IfcSpace')
        in_building=[ space.id() not in list(self._building_elements.keys()) for space in ifcspaces]
        
        space_shapes = [ifcelement_as_solid(x) for x in compress(ifcspaces,in_building) if x.Representation is not None]
        #print(ifcspaces)
        #print(space_shapes)
        #print([self._building]+space_shapes)
        
        solids=shapes_as_solids([self._building]+space_shapes)
        #print(solids)
        #filled_building=fuse_listOfShape(solids)
        
        self._external_shell = get_external_shell(solids)
        
        
        wall_shapes={ el.id():building_shapes[el.id()] for el in self._building_elements.values() if el.is_a()=='IfcWall'}
        self._wall_norm_plane={}
        for id,wall_shape in wall_shapes.items():
            self._wall_norm_plane[id]=exterior_wall_normal_and_plane(wall_shape,self._external_shell)
        
    def display(self,disp_external_shell=False,disp_solar_shapes=False):
        
        
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
        
        
        self._results=dict()
        
        for id,facelist in self._solar_faces.items():
            
            #if id!=836:
            #    continue
            idhost = self._hosting_solar_id[id] 
            norm,plane =self._wall_norm_plane[idhost]
            
            sun_to_earth_project,irradiance=self._proj_loc.data_critical_period_day(facelist[0])
            
            print(" Opening ID  : ", id)
            rtaa=rtaa_on_faces(facelist,plane,self._building,sun_to_earth_project,cut,adjust)
            #rtaa.compute_masks_hangover(hp=.85,lp=.5,dhp=.15,dhm=x[k])
            
            rtaa.compute_masks()
            rtaa.compute_cm(irradiance)
            self._results[id]=rtaa
            #rtaa.display()
            
        
    def cm(self):
        return { id : r._cm for id,r in self._results.items()} 
        
    def raw_result(self):
        #return dict of dataframes
        pass

class rtaa_ventilation_study:
    
    def __init__(self,ifcfilename):
        setting=ifcopenshell.geom.settings()
        setting.set(setting.USE_PYTHON_OPENCASCADE, True)
        self._ifc_file= ifcopenshell.open(ifcfilename)
        
        self._space_elements=dict()
        self._opening_elements=dict()
        
        print(" ifc file ",self._ifc_file)
        
        self._proj_loc=project_location()
        self._proj_loc.set_location_from_ifc(self._ifc_file)
        self._proj_loc.set_northing_from_ifc(self._ifc_file)

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
          
        
        if container=="spaces":
            self._space_elements={el.id():el for el in element_set}
            
        elif container=="openings":
            self._opening_elements={el.id():el for el in element_set}
        
    
    
    def _remove_elements_by_ids(self,ids=[],container=None):
        if(len(ids)==0):
            print(" Error : provides ids or classes to set elements of interest")
            return
        if container is None:
            print(" Error : provide a container")
            
                    
        if container=="spaces":
            for id in ids:
                self._space_elements.pop(id)
        elif container=="openings":
            for id in ids:
                self._opening_elements.pop(id)
    
    
    def add_space_elements(self,ids=[],types=[]):
        self._add_elements(ids,types,'spaces')
        print(" Number of elements for solar analysis: ",len(self._space_elements.keys()))
    
    def remove_spaces_elements(self,ids=[]):
        self._remove_elements(ids,types,'spaces')
        
    
    
    
        
    def add_opening_elements(self,ids=[],types=[]):
        self._add_elements(ids,types,'openings')
        print(" Number of elements in building: ",len(self._opening_elements.keys()))
    
    def remove_openings_elements(self,ids=[]):
        self._remove_elements(ids,types,'openings')    
                

    def set_geometries(self):
        
        op_ids=set(self._opening_elements.keys())
        # for opening of interest only 
        # a subset of wall geometry and windows
        self._window_by_wall={}
        self._wall_shapes={}
        self._window_shapes={}
        
        ifcwalls= self._ifc_file.by_type('IfcWall')
        
        for w in ifcwalls:
            l_op= opening_in_wall(w)
            op_in_wall= set(l_op)
            if len(op_in_wall.intersection(op_ids))>0:
                
                self._window_by_wall[w.id()]=l_op
                self._wall_shapes[w.id()]=ifcelement_as_solid(w)
                
                for op in l_op:
                    opening=self._ifc_file.by_id(op)
                    self._window_shapes[opening.id()]=ifcelement_as_solid(opening)
    
    def display(self):
        
        
        def rgb_color(r, g, b):
            return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
        x=50/256
        gray=rgb_color(x, x, x)
        orange=rgb_color(255/256,127/256,80/256)
        display, start_display, add_menu, add_function_to_menu = init_display()
        
        ifcspaces=[self._ifc_file.by_id(id) for id in self._space_elements]
        spaceshapes=[ ifcelement_as_solid(s) for s in ifcspaces]
        
        [display.DisplayShape(s,color=orange,transparency=0.1)for s in spaceshapes]
        [display.DisplayShape(s,color=gray,transparency=0.8)for s in self._wall_shapes.values()]
        
        gpp=GProp_GProps()
        
        for svd in self._results:
            
            #print(svd._wall_faces)
            for f in svd._wall_faces.values():
                display.DisplayShape(f,color='BLUE',transparency=0.1)
                
            for f in svd._win_faces.values():
                display.DisplayShape(f,color='GREEN',transparency=0.1)    
                
        
        display.FitAll()
        
        start_display()
        
    def export_raw(self):
        print([hasattr(r,'_win_data') for r in self._results])
        win=pd.concat([r._win_data for r in self._results if hasattr(r,'_win_data')])
        wall = pd.concat([r._wall_data for r in self._results if hasattr(r,'_wall_data')])
        
        return win,wall
           
    
    def run(self):
        self._results=[]
        for ifcspaceid in self._space_elements:
            
            ifcspace = self._ifc_file.by_id(ifcspaceid)
            svd=space_ventilation_data(ifcspace)
            #print("  Processing ",ifcspace)
            
            svd.extract_faces(ifcspace,
                            self._window_by_wall,
                            self._wall_shapes,
                            self._window_shapes)
            #print("\n\n")
            svd.info()
            svd.sweeping(self._window_shapes)
            #svd.analyze_faces()
            svd.opening_ratio(self._proj_loc)
            
            self._results.append(svd)
        win,wall=self.export_raw()
        print(win)
        print(wall)
       
            

class space_ventilation_data:
    def __init__(self,ifcspace):
        self._space=ifcspace
        self._win_by_wall=defaultdict(list)
        self._wall_faces =defaultdict(list)
        self._win_faces  =defaultdict(list)
        

    
    def extract_faces(self,ifcspace,window_by_wall,wall_shapes,windows_shapes):
        if ifcspace.Representation is not None:
            ss =create_shape(setting, ifcspace).geometry
        else :
            print(" No geometry for ifcspace : ",ifcspace)
            return
        
        
        shu=ShapeUpgrade_UnifySameDomain(ss)
        shu.Build()
        ss=shu.Shape()
        
        origin = gp_Pnt(0.0,0.0,0.0)
        Zaxis = gp_Ax1(origin,gp_Dir(0.0,0.0,1.0))
        
        faces=list(TopologyExplorer(ss).faces())
        # need to fuse face that have commom plane
        lax3=[]
        lface=[]
        lnorm=[]
        for f in faces:
            
            srf = BRep_Tool().Surface(f)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(f.Orientation()==1):
                face_norm.Reverse()
        
            lax3.append(plane.Position())
            lface.append(f)
            lnorm.append(face_norm)
            #print('extrusion ',ext_vec.Coord())
            
            
            # skipping horizontal face but referencing the soil one
            if( face_norm.Dot(Zaxis.Direction())<(-1+1e-5)):
                self._soil_face=f
                continue
        #make this test at the end to enlarge the face ?
        lcoplanar=[[] for x in lax3 ]
        for  i,j in itertools.combinations(range(0,len(lax3)),r=2):
            coplanar= lax3[i].IsCoplanar(lax3[j],1.e-5,1.e-5)
            if coplanar:
                lcoplanar[i].append(j)
                print(i,' ',j,' ')#,lcoplanar[i])
        newfaces=[]
        toremove=[]
        for f,indices in zip(faces,lcoplanar):
            if len(indices)>0:
                los=[f]
                [los.append(lface[i]) for i in indices]
                toremove.extend(indices)
                print(los)
                new=fuse_listOfShape(los)
                #new=list(TopologyExplorer(new).faces())
                print('new ',new)
                newfaces.append(new)
            else:
                newfaces.append(f)
        # remove face that are in indices 
        toremove.sort(reverse=True)
        #print(' new face ',len(newfaces))
        #print(' to remove ',toremove)
        for i in toremove:
            newfaces.pop(i)
            lnorm.pop(i)
            lcoplanar.pop(i)
                
        #print(' new face ',newfaces)
        for f,face_norm,coplanar in zip(newfaces,lnorm,lcoplanar):    
            ext_vec=gp_Vec(face_norm)*.1
            extrusion = BRepPrimAPI_MakePrism(f,ext_vec,False,True).Shape()
            
            #linter.append(extrusion)
            
            for wall_id,lwin in window_by_wall.items():
                wall=wall_shapes[wall_id]
                
                intersection=BRepAlgoAPI_Common(extrusion,wall)
                intersection_wall=intersection.Shape()
                intersection_wall_solids=list(TopologyExplorer(intersection_wall).solids())
                            
                #the face of space catch a wall shape
                if len(intersection_wall_solids)>0:
                    # Searching if a window is hosted by the 
                    # portion of wall catched by face extrusion                                         
                    for win_id in lwin:
                        #print('win_id ',win_id)
                        win=windows_shapes[win_id]
                        
                        intersection=BRepAlgoAPI_Common(extrusion,win)
                        intersection_win=intersection.Shape()
                        intersection_win_solids=list(TopologyExplorer(intersection_win).solids())
                        # the wall face catch a window
                        # extracting datas from window
                        if len(intersection_win_solids)>0:
                            bigfaces=biggestface_along_vector(win,face_norm)
                            #lbigface.extend(bigfaces)
                            #lface.append(f)
                            #faceswin[win_id].extend(bigfaces)
                            self._win_by_wall[wall_id].append(win_id)
                            #self._wall_faces[wall_id].append(f)
                            self._wall_faces[wall_id].extend(TopologyExplorer(f).faces())
                            print(wall_id,' ',self._wall_faces[wall_id])
                            self._win_faces[win_id].extend(bigfaces)
                            #svd.update(wall_id,win_id,f,bigfaces)
    
    
    def opening_ratio(self,projloc):
        if(len(self._win_by_wall.keys())==0):
            print("No opening in this space")
            return
        #inverse mapping
        wall_by_win={}
        for wall_id,win_ids in self._win_by_wall.items():
            for win_id in win_ids:
                wall_by_win[win_id]=wall_id
        
        gpp=GProp_GProps()
        
        
        wall_area=defaultdict(list)
        wall_area_total=dict()
        wall_angle=dict()
        wall_orient=dict()
        for wall_id,f_list in self._wall_faces.items():
            angles=[]
            orient=[]
            for f in f_list:
                brepgprop_SurfaceProperties(f,gpp)
                wall_area[wall_id].append(gpp.Mass())
                angles.append(round(projloc.face_orientation_angle_tn(f),6))
                orient.append(projloc.face_orientation_sector(f))
            wall_angle[wall_id]=list(set(angles))[0]
            wall_orient[wall_id]=list(set(orient))[0]
            wall_area_total[wall_id]=sum(wall_area[wall_id])
        
        
        win_area=defaultdict(list)
        win_area_total=dict()
        win_area_by_wall=defaultdict(float)
        
        
        win_rtaa_poro=dict()
                       
        for win_id,f_list in self._win_faces.items():
            win_rtaa_poro[win_id]=1.0
            for f in f_list: 
                brepgprop_SurfaceProperties(f,gpp)
                
                win_area[win_id].append(gpp.Mass())
            win_area_total[win_id]= sum(win_area[win_id]) 
            
            win_area_by_wall[wall_by_win[win_id]]+=sum(win_area[win_id])
        
        # room winid poro
        win_ids=list(self._win_faces.keys())
        space_name= [self._space.Name]*len(win_ids)
        space_id= [self._space.id()]*len(win_ids)
        poros= [ win_rtaa_poro[win_id] for win_id in  win_ids]
        areas =[ win_area_total[win_id] for win_id in  win_ids]
        
        data={'space_name':space_name,
              'space_id':space_id,
              'win_ifc_id':win_ids,
              'rtaa_poro':poros,
              'area':areas}
        win_data=pd.DataFrame(data)
        self._win_data=win_data
        
        print(win_data)
        
        wall_ids=list(self._win_by_wall.keys())
        space_name= [self._space.Name]*len(wall_ids)
        space_id= [self._space.id()]*len(wall_ids)
        wall_areas =[ wall_area_total[wall_id] for wall_id in  wall_ids]
        win_areas =[ win_area_by_wall[wall_id] for wall_id in  wall_ids]
        wall_angle=[wall_angle[wall_id] for wall_id in wall_ids]
        wall_orient=[ wall_orient[wall_id] for wall_id in wall_ids]
        data={'space_name':space_name,
              'space_id':space_id,
              'wall_ifc_id':wall_ids,
              'wall_area':wall_areas,
              'win_area':win_areas,
              'wall_angle':wall_angle,
              'wall_orient':wall_orient}
        wall_data=pd.DataFrame(data)
        wall_data.sort_values(by=['win_area'],ascending=False,inplace=True)
        self._wall_data=wall_data
        print(wall_data)
        
        self._opening_ratio=wall_data['win_area'].sum()/wall_data['wall_area'][0]
        print(self._opening_ratio)
        
        
        """
        largest_opened_wall_id=max(win_area_by_wall,key=win_area_by_wall.get)
        # largest window area
        print('\n\n')
        print('Area of each window ')
        for k,v in win_area_total.items():
            print('  win id ',k,' : ',v)
        #print(win_area_total)
        print(' Area of window on each wall ')
        for k,v in win_area_by_wall.items():
            print('  wall id ',k,' : ',v, ' // wall area : ',wall_area_total[k] )
       
        print(" largest opening is ",largest_opened_wall_id)
        #for k,v in self._win_by_wall.items():
        #    if largest_windows_id in v:
        #        wall_largest_window=k
        
        # print the details : Area of each opening, associated wall area
        
        other_walls=list(self._win_by_wall.keys())
        other_walls.remove(largest_opened_wall_id)
        other_windows_area=list()
        for wall_id in other_walls:
            for win_id in self._win_by_wall[wall_id]:
                # aera of the window / area of the hosting wall
                ratio = win_area_total[win_id]/wall_area_total[ wall_by_win[win_id]]
                #print( win_id,' ',ratio)
                if ratio >0.1:
                    other_windows_area.append(win_area_total[win_id])
        result={}
        result['A1']= wall_area_total[largest_opened_wall_id]
        result['A2']= win_area_by_wall[largest_opened_wall_id]
        result['A3s']= other_windows_area
        print(' Opening Analysis results ')
        print(result)
        print('     opening ratio ', (result['A2']+sum(result['A3s']))/result['A1'])
      
        
        porosity=1.0
        """
    def sweeping(self,windows_shapes):
        # compute and analyze intersection face and curve linking two windows
        
        if( len(self._win_by_wall.keys())==0):
            print("No window in this space")
            return
        
        if( len(self._win_by_wall.keys())==1):
            print("All windows hosted by the same wall , no need to compute sweeping")
            return
        
        # preprocessing of the soil surface
        gpp=GProp_GProps()
        srf = BRep_Tool().Surface(self._soil_face)
        plane = Geom_Plane.DownCast(srf)
                    
        # convert edges to lines with min max value of paramter
        edges=list(TopologyExplorer(self._soil_face).edges())
        llines_face=[]
        for e in edges:
            (curve,minu,maxu)=BRep_Tool.Curve(e)
            adapt=GeomAdaptor_Curve(curve)
            llines_face.append( (projlib_Project(plane.Pln(),adapt.Line()),minu,maxu))
        
        # convert vertices to points            
        vertices=list(TopologyExplorer(self._soil_face).vertices())
        soil_pt=[BRep_Tool.Pnt(v) for v in vertices]
        
        # computation of an approximated diagonal
        brepgprop_SurfaceProperties(self._soil_face,gpp)
        mc_soil=gpp.CentreOfMass()
        #distance from the mass center is already a proxy for half diag
        half_diag = sum([mc_soil.Distance(p) for p in soil_pt])/len(soil_pt)
        
        #print('Demie diagonale ', half_diag)
        
        # mass center of each windows for this room
        mass_centers={}
        for wall_id,win_ids in self._win_by_wall.items():
            srf_wall=BRep_Tool().Surface(self._wall_faces[wall_id][0])
            
            for win_id in win_ids:
                brepgprop_VolumeProperties(windows_shapes[win_id],gpp)
                mc=gpp.CentreOfMass()
                proj=GeomAPI_ProjectPointOnSurf(gpp.CentreOfMass(),srf_wall)
                mass_centers[win_id]=proj.Point(1)
        
        length_between_mc={}
        for win_id1,win_id2 in combinations(mass_centers.keys(),2):
            mc1=mass_centers[win_id1]
            mc2=mass_centers[win_id2]
            #print(" Distance between ",win_id1,' ',win_id2)
        #for i,mc1 in enumerate(mass_centers[:-1]):
        #    for mc2 in mass_centers[i+1::]:
                #print(" new line ")
                #print( mc.Coord(),' ',mc2.Coord())
                
            #line2d between the two mass centers
            mc12d = projlib_Project(plane.Pln(),mc1)
            mc22d = projlib_Project(plane.Pln(),mc2)
            between_vec = gp_Vec2d(mc22d,mc12d).Normalized()
            lin2d=gp_Lin2d( mc12d, gp_Dir2d(between_vec))
            
            # intersection between the line and the soil_face edges (as line)                    
            lin2d_values=[]
            for l,minu,maxu in llines_face:
                inter=Geom2dAPI_InterCurveCurve(Geom2d_Line(l),Geom2d_Line(lin2d),1e-5)
                #print(" number of intersections ",inter.NbPoints())
                if inter.NbPoints()>0:
                    uvalue=elclib_Parameter(l,inter.Point(1))
                    if( (uvalue>=minu) & (uvalue<=maxu)):
                        uvalue2=elclib_Parameter(lin2d,inter.Point(1))
                        lin2d_values.append(uvalue2)
                #else :
                    # rare case of line joining two window parallel to another wall
                    
                
                
                # take only the intersection inside the surface

            
            # more intersection than just the two window mc
            intersection_points=[]
            if len(lin2d_values)>2:
                #print(lin2d_values)
                                        
                lmax= elclib_Parameter(lin2d,mc12d)
                lmin= elclib_Parameter(lin2d,mc22d)
                #print(lmin,' ',lmax)
                
                #values_inter=[ v for v in lin2d_values if (v>lmax) | (v<lmin)]
                
                for v in lin2d_values:
                    #print(v) 
                    pt2d= elclib_Value(v,lin2d)
                    
                    # lower bound of the line                            
                    if (pt2d.IsEqual(mc12d,1e-5)) :
                        continue
                    # higher bound of the line    
                    if (pt2d.IsEqual(mc22d,1e-5)):
                        continue
                    
                    if( (v >= lmax) | (v<=lmin)):
                        continue
                    #print("good v ") 
                    pt3d= elclib_To3d(plane.Pln().Position().Ax2(),pt2d)
                    intersection_points.append(pt3d)
                    
        
            # Only one intersection mean that the ray do not get back inside
            if len(intersection_points)==1:
                continue
            
                            
            #no problem if no intersection 
            to_keep_idx=set()
            for p in intersection_points:
                d=[v.Distance(p) for v in soil_pt]
                to_keep_idx.add( d.index(min(d)))
            to_keep=[soil_pt[i] for i in to_keep_idx]    
            
            path=[mc1]+to_keep+[mc2]
            length=0.0
            for p1,p2 in pairwise(path):
                length+=p1.Distance(p2)
                #print('lenth ',length ,' divided (must be >1) ', length/(half_diag))
                length_between_mc[(win_id1,win_id2)]=length/half_diag
        
        print(' Sweeping (must be >1) : ' )
        for k,v in length_between_mc.items():
            print('     ',k[0],' ',k[1],' ',v)
        print(length_between_mc)
        return length_between_mc
    
      
    def info(self):
        print('Space Name ',self._space.Name)
        print('Space Id   ',self._space.id())
        print('*** Windows id by wall id of this space')
        for wall_id,win_list in self._win_by_wall.items():
            print("     wall id ",wall_id)
            for w in win_list:
                print("         win id ",w)
        print('*** Space face by wall id in this space')
        for wall_id,face_list in self._wall_faces.items():
            print("     wall id ",wall_id)
            for w in face_list:
                print("         face ",w)
        print('*** Largest faces of window in this space')
        for wall_id,face_list in self._win_faces.items():
            print("     wall id ",wall_id)
            for w in face_list:
                print("         face ",w)
        print('\n')            
     

if __name__ == "__main__":

    # due to some bugs in ipython parsing
    __import__("logging").getLogger("parso.python.diff").setLevel("INFO")
    __import__("logging").getLogger("parso.cache").setLevel("INFO")
    __import__("logging").getLogger("asyncio").setLevel("INFO")

    

    # Initialize a graphical display window (from ifcos)

    

    filename = 'data/Rtaa_validation_run_onewindow.ifc'
    filename = 'data/Rtaa_validation_run.ifc'
    filename = 'data/Rtaa_validation_run_joues2.ifc'
    #ifc_file= ifcopenshell.open('data/simple_reunion_northaligned.ifc')
    #ifc_file= ifcopenshell.open('data/Rtaa_validation_run.ifc')


    filename='C:/Users/cvoivret/source/ventilation_test/DCE_CDV_BAT.ifc'
    
    #filename='DCE_CDV_BAT.ifc'

    rsv=rtaa_ventilation_study(filename)
    lgt4b=[499,705,731,757,783,809,859]
    windows=[14407,15736,7565]
    
    rsv.add_space_elements([],['IfcSpace'])#,['IfcSpace'])
    rsv.add_opening_elements([],['IfcWindow'])
    rsv.set_geometries()
    rsv.run()
    rsv.display()
    
    
    
    
    """
    tn_angle=[0.0,np.pi*.5,np.pi,3.*np.pi*.5]
    orientations_name=['NORD','EST','SUD','OUEST']
    
    config={name:(angle) for name,angle in zip(orientations_name,tn_angle)}
        
    res=[]
    
    totake=orientations_name#['NORD']
    list_angle= [ config[n] for n in totake]
    
    x=np.arange(0.15,1.15,.1)
    
    filename = 'data/debords_casquettes_fins_triangle.ifc'
    rss=rtaa_solar_study(filename)
    rss.add_building_elements([],['IfcWall','IfcSlab'])
    rss.add_solar_elements([],['IfcWindow','IfcDoor'])
    rss.set_geometries()
        
    cm_orientation=defaultdict(list)   
    for angle,orient in zip(list_angle,totake):
        rss._proj_loc.update_northing_from_angle(angle)
        #rss.display(True)
        rss.run()
        cm_orientation[orient]=rss.cm()
    
    res=pd.DataFrame(cm_orientation)
    res=res.sort_values(by='NORD',ascending=False)
    res.to_csv(filename[:-3]+'csv')
    
    """
    
    """
    filename_list=['data/debords_casquettes_epais.ifc',
                    'data/debords_casquettes_fins.ifc',
                    'data/debords_epais.ifc',
                    'data/debords_fins.ifc']
                    
                    
                    
    #filename = 'data/debords_casquettes_epais.ifc'
    
    for filename in filename_list:
        rss=rtaa_solar_study(filename)
        rss.add_building_elements([],['IfcWall','IfcSlab'])
        rss.add_solar_elements([],['IfcWindow','IfcDoor'])
        rss.set_geometries()
            
        cm_orientation=defaultdict(list)   
        for angle,orient in zip(list_angle,totake):
            rss._proj_loc.update_northing_from_angle(angle)
            #rss.display(True)
            rss.run()
            cm_orientation[orient]=rss.cm()
        
        res=pd.DataFrame(cm_orientation)
        res=res.sort_values(by='NORD',ascending=False)
        res.to_csv(filename[:-3]+'csv')
    
    rtaa_data= pd.read_excel('data/debords_casquette.xlsx')
        
    colors=['r','g','b','k']
    marks=['x','d','o']
    #for cmo,m in zip(cm_orientation,marks):
    for name,c in zip(cm_orientation.keys(),colors):
        data=list(cm_orientation[name].values())
        data.sort(reverse=True)
        plt.plot(x,data,'-',color=c,marker='o',lw=2,label=name+'_computed')
    
    for name,c in zip(totake,colors):
        plt.plot(x,rtaa_data[name],'--',color=c,lw=2,label=name+'_ref')
    plt.ylim(0,1)
    plt.xlim(0.1,1.1)
    plt.legend()
    plt.show()
    
    """
   
    #filename = 'data/villa.ifc'

    #ifc_file= ifcopenshell.open(filename)
    
    """
    rsv=rtaa_ventilation_study(filename)
    lgt4b=[499,705,731,757,783,809,859]
    lgt4b=[499]

    
    rsv.add_space_elements(lgt4b,[])#,['IfcSpace'])
    rsv.add_opening_elements([],['IfcWindow','IfcDoor'])
    
    
    #rsv.add_space_elements([],['IfcSpace'])
    #rsv.add_opening_elements([],['IfcWindow'])
    
    rsv.set_geometries()
    rsv.run()
    rsv.display()
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
    
    
