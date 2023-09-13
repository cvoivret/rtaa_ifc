from collections import defaultdict
from itertools import tee,combinations
import ifcopenshell
from ifcopenshell.geom import create_shape
from ifcopenshell.geom.occ_utils import yield_subshapes

from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color,Quantity_TOC_RGB


from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox,BRepPrimAPI_MakePrism,BRepPrimAPI_MakeHalfSpace,BRepPrimAPI_MakeSphere,BRepPrimAPI_MakeCylinder
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties,brepgprop_VolumeProperties,brepgprop_LinearProperties
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing,BRepBuilderAPI_MakeSolid,BRepBuilderAPI_MakeFace,BRepBuilderAPI_MakeEdge
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

from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.Precision import precision


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
from OCC.Core.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from OCC.Core.ProjLib import ProjLib_ProjectOnPlane
from OCC.Core.ProjLib import projlib_Project

from OCC.Core.GeomProjLib import geomprojlib_Project
from OCC.Core.GeomProjLib import geomprojlib_ProjectOnPlane

from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve

from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector

from OCC.Core.Standard import standard_Purge

from OCC.Extend.DataExchange import write_stl_file

from OCC.Core.AIS import AIS_Line
from OCC.Core.GeomAPI import geomapi_To3d
from OCC.Core.Geom2d import      Geom2d_Line
from OCC.Core.ElCLib import elclib_Parameter,elclib_To3d,elclib_Value

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

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

def ifcelement_as_solid(ifcwin):
        if ifcwin.Representation is None:
            return None
        shape=create_shape(setting, ifcwin).geometry
        solids=shapes_as_solids([shape])
        solid=fuse_listOfShape(solids)
        return solid

class space_ventilation_data:
        def __init__(self,ifcspace):
            self._space=ifcspace
            self._win_by_wall=defaultdict(list)
            self._wall_faces =defaultdict(list)
            self._win_faces  =defaultdict(list)

        #def update(self,wall_id,win_id,wall_face,win_face):
        #    self._win_by_wall[wall_id].append(win_id)
        #    self._wall_faces[wall_id].append(wall_face)
        #    self._win_faces[win_id].extend(win_face)
            
        
        
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
            for f in faces:
                soil_face=False
                srf = BRep_Tool().Surface(f)
                plane = Geom_Plane.DownCast(srf)
                face_norm = plane.Axis().Direction()
                if(f.Orientation()==1):
                    face_norm.Reverse()
                ext_vec=gp_Vec(face_norm)*.1
                #print('extrusion ',ext_vec.Coord())
                
                # skipping horizontal face but referencing the soil one
                if( face_norm.Dot(Zaxis.Direction())<(-1+1e-5)):
                    self._soil_face=f
                    continue
                
                
                extrusion = BRepPrimAPI_MakePrism(f,ext_vec,False,True).Shape()
                #linter.append(extrusion)
                #print('   ')
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
                                self._wall_faces[wall_id].append(f)
                                self._win_faces[win_id].extend(bigfaces)
                                #svd.update(wall_id,win_id,f,bigfaces)
      
            
        def opening_ratio(self):
            if(len(self._win_by_wall.keys())==0):
                print("No windows in this space")
                return
            #inverse mapping
            wall_by_win={}
            for wall_id,win_ids in self._win_by_wall.items():
                for win_id in win_ids:
                    wall_by_win[win_id]=wall_id
            
            gpp=GProp_GProps()
            
            
            wall_area=defaultdict(list)
            wall_area_total=dict()
            win_by_wall_area_total=dict()
            for wall_id,f_list in self._wall_faces.items():
                for f in f_list:
                    brepgprop_SurfaceProperties(f,gpp)
                    wall_area[wall_id].append(gpp.Mass())
                wall_area_total[wall_id]=sum(wall_area[wall_id])
            
            win_area=defaultdict(list)
            win_area_total=dict()
            win_area_by_wall=defaultdict(float)
            for win_id,f_list in self._win_faces.items():
                for f in f_list: 
                    brepgprop_SurfaceProperties(f,gpp)
                    win_area[win_id].append(gpp.Mass())
                win_area_total[win_id]= sum(win_area[win_id]) 
                
                win_area_by_wall[wall_by_win[win_id]]+=sum(win_area[win_id])
                
            # largest window area
            print(win_area_total)
            largest_opened_wall_id=max(win_area_by_wall)
            #for k,v in self._win_by_wall.items():
            #    if largest_windows_id in v:
            #        wall_largest_window=k
            
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
            print(result)
            print(' opening ratio ', (result['A2']+sum(result['A3s']))/result['A1'])
            """
            opening_ratio=(win_area_total[largest_windows_id]+other_windows_area)/wall_area_total[wall_largest_window]
            print(" wall area of the largest window",wall_largest_window," ",wall_area_total[wall_largest_window])
            print(" largest window area",win_area_total[largest_windows_id])
            print(" sum of others window area ", other_windows_area)
            print(" opening ratio (no porosity) ",opening_ratio)
            """
            
            porosity=1.0
        
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
            
            print('Demie diagonale ', half_diag)
            
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
                print(" Distance between ",win_id1,' ',win_id2)
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
                    inter=Geom2dAPI_InterCurveCurve(Geom2d_Line(l),Geom2d_Line(lin2d))
                    uvalue=elclib_Parameter(l,inter.Point(1))
                    # take only the intersection inside the surface
                    if( (uvalue>=minu) & (uvalue<=maxu)):
                        uvalue2=elclib_Parameter(lin2d,inter.Point(1))
                        lin2d_values.append(uvalue2)
                
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
                    print('lenth ',length ,' divided (must be >1) ', length/(half_diag))
                    length_between_mc[(win_id1,win_id2)]=length
                        
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
            #print(' Largest faces of window in this space\n',self._win_faces)
                

if __name__ == "__main__":

    # due to some bugs in ipython parsing
    __import__("logging").getLogger("parso.python.diff").setLevel("INFO")
    __import__("logging").getLogger("parso.cache").setLevel("INFO")
    __import__("logging").getLogger("asyncio").setLevel("INFO")

        
    setting=ifcopenshell.geom.settings()
    setting.set(setting.USE_PYTHON_OPENCASCADE, True)

    filename = 'data/villa.ifc'

    ifc_file= ifcopenshell.open(filename)
    
    ifcwalls=ifc_file.by_type('IfcWall')
    ifcspaces=ifc_file.by_type('IfcSpace')
    ifcwindows=ifc_file.by_type('IfcWindow')

    # linking wallid to a list of window id
    window_by_wall1={}
    for w in ifcwalls:
        l= window_in_wall(w)
        if len(l)>0:
            window_by_wall1[w.id()]=l
    # collecting ifcwall object that host a window
    ifc_wall_with_windows = [ ifc_file.by_id(id) for id in window_by_wall1.keys()]
     
    # creating shape and store by object if in a dict
    windows_shapes1 = { x.id():ifcelement_as_solid(x) for x in ifcwindows if x.Representation is not None }
    wall_shapes1  = {x.id():   create_shape(setting, x).geometry for x in ifc_wall_with_windows if x.Representation is not None}
    
    for ifcspace in ifcspaces:
        svd=space_ventilation_data(ifcspace)
        svd.extract_faces(ifcspace,window_by_wall1,wall_shapes1,windows_shapes1)
        print("\n\n")
        svd.info()
        svd.sweeping(windows_shapes1)
        svd.opening_ratio()
    
   
    
    
    #Îwall_shapes2  = [create_shape(setting, x).geometry for x in ifcwalls if x.Representation is not None]
    #space_shapes = [create_shape(setting, x).geometry for x in ifcspaces if x.Representation is not None]
    
    #lsphere=[  	BRepPrimAPI_MakeSphere(p, 1.0).Shape() for p in lpoints]

    
    """
    def rgb_color(r, g, b):
        return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    
    x=50/256
    gray=rgb_color(x, x, x)

    display, start_display, add_menu, add_function_to_menu = init_display()
    context=display.GetContext()

    
    #display.DisplayShape(wall_shapes,color=gray,transparency=0.9)
    #[display.DisplayShape(s,color='BLUE',transparency=0.1) for s in linter]
    #[display.DisplayShape(s,color='RED',transparency=0.9) for s in lface]
    [display.DisplayShape(s,color='GREEN',transparency=0.1) for s in lbigface]

    #[display.DisplayShape(s,color=gray,transparency=0.9) for s in wall_shapes2]
    [display.DisplayShape(s,color='BLUE',transparency=0.9) for s in space_shapes]
    #[context.Display(AIS_Line(s),True) for s in llines]
    [display.DisplayShape(BRepPrimAPI_MakeSphere(c,0.1).Shape() )for c in lpoints]
    display.FitAll()
    start_display()
    """
    
    