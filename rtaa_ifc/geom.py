
import ifcopenshell
from ifcopenshell.geom import create_shape
from ifcopenshell.geom.occ_utils import yield_subshapes

setting=ifcopenshell.geom.settings()
setting.set(setting.USE_PYTHON_OPENCASCADE, True)

from collections import defaultdict
from collections import Counter

from OCC.Core.BRep import BRep_Tool

from OCC.Core.TopoDS import TopoDS_Face,TopoDS_Shape
from OCC.Core.TopTools import TopTools_ListOfShape,TopTools_IndexedMapOfShape
from OCC.Core.TopExp import topexp,TopExp_Explorer 
from OCC.Core.ShapeExtend import ShapeExtend_Explorer
from OCC.Core.TopAbs import TopAbs_SOLID,TopAbs_FACE,TopAbs_SHELL,TopAbs_WIRE,TopAbs_SHAPE 
import OCC.Core.ShapeFix as ShapeFix_Shape
from OCC.Core.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCC.Core.BOPAlgo import BOPAlgo_BOP,BOPAlgo_Operation
from OCC.Core.BOPTools import BOPTools_AlgoTools

from OCC.Core.GProp import GProp_GProps
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.BRepGProp import brepgprop #SurfaceProperties,VolumeProperties
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing,BRepBuilderAPI_MakeSolid,BRepBuilderAPI_MakeFace
#from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Copy,BRepBuilderAPI_Transform
from OCC.Core.Geom import Geom_Plane

#from OCC.Core.BOPAlgo import BOPAlgo_CellsBuilder

from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox#,BRepPrimAPI_MakePrism,BRepPrimAPI_MakeHalfSpace,BRepPrimAPI_MakeSphere,BRepPrimAPI_MakeCylinder
from OCC.Extend.TopologyUtils import TopologyExplorer#, WireExplorer,TopoDS_Iterator


def face_mesh_triangle(comp=TopoDS_Shape(), isR=0.1, thA=0.1):
    # Mesh the shape
    BRepMesh_IncrementalMesh(comp, isR, True, thA, True)
    bild1 = BRep_Builder()
    comp1 = TopoDS_Compound()
    bild1.MakeCompound(comp1)
    bt = BRep_Tool()
    ex = TopExp_Explorer(comp, TopAbs_FACE)
    while ex.More():
        face = topods_Face(ex.Current())
        location = TopLoc_Location()
        facing = bt.Triangulation(face, location)
        tab = facing.Nodes()
        tri = facing.Triangles()
        print(facing.NbTriangles(), facing.NbNodes())
        for i in range(1, facing.NbTriangles() + 1):
            trian = tri.Value(i)
            index1, index2, index3 = trian.Get()
            for j in range(1, 4):
                if j == 1:
                    m = index1
                    n = index2
                elif j == 2:
                    n = index3
                elif j == 3:
                    m = index2
                me = BRepBuilderAPI_MakeEdge(tab.Value(m), tab.Value(n))
                if me.IsDone():
                    bild1.Add(comp1, me.Edge())
        ex.Next()
    return comp1
    
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

# try to convert a shape in a list of solid 
# if needed decomposing a compound
# make a solid from a list of faces
# handle compound (hopefully with solid) and a faces as shell

def solids_from_shape(shape):
    lsolid=[]
    print(' solids_from_shape ----',shape) 
    maps=TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape,TopAbs_SOLID,maps)
    # extract solids
    maps.Clear()
    if(maps.Size()>0):
        #print('number of solids ', maps.Size())
        lsolid.extend([maps.FindKey(i) for i in range(1,maps.Size()+1)])
    # extract faces ans sew them
    else:
        maps.Clear()
        topexp.MapShapes(shape,TopAbs_FACE,maps)
        if( maps.Size()>0):
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


def shapes_as_solids(lshape):
    """
    Try to build a list of solid from a list to shapes.
    Flatten nested shapes.
    Sew tesselated shape to build a solid.
    Fix shape if needed.
    
    """
    lsolid=[]
    print(' shape as solid ----',lshape) 
    maps=TopTools_IndexedMapOfShape()
    for s in lshape:
        #print(' shape ',s)
        maps.Clear()
        topexp.MapShapes(s,TopAbs_SOLID,maps)
        se=ShapeExtend_Explorer()
        lsh=se.SeqFromCompound(s,True)
        #print(lsh)
        
        #print('list of shapes --- ')
        #for i in  range(1,lsh.Length()+1):
        #    print('shape ', lsh.Value(i))
        
        if(maps.Size()>0):
            #print('number of solids ', maps.Size())
            lsolid.extend([maps.FindKey(i) for i in range(1,maps.Size()+1)])
        else:
            maps.Clear()
            topexp.MapShapes(s,TopAbs_FACE,maps)
            if( maps.Size()==0):
                continue
            #print(maps.Size())
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
def ifcelement_as_solid_old(element):
    if element.Representation is None:
        return None
    shape=create_shape(setting, element).geometry
    print('\n ifc shape ',shape)
    solids=shapes_as_solids([shape])
    print('list of solids ',solids)
    solid=None
    if( len(solids)==1):
        solid = solids[0]
    else :
        solid = fuse_listOfShape(solids)
    print('solid after fuse: ',solid)
    if (solid is not None ) :
        shu=ShapeUpgrade_UnifySameDomain(solid)
        shu.Build()
        solid=shu.Shape()    
    #print(solid)
    
    # add unifiysamedomain ?
    
    return solid
    
def ifcelement_as_solid(element):
    if element.Representation is None:
        return None
        
    shape=create_shape(setting, element).geometry
    
    solids=solids_from_shape(shape)
    solid=None
    if( len(solids)==1):
        solid = solids[0]
    else :
        solid = fuse_listOfShape(solids)
        
    print('solid after fuse: ',solid)
    
    if (solid is not None ) :
        shu=ShapeUpgrade_UnifySameDomain(solid)
        shu.Build()
        solid=shu.Shape()    
    
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
    
    BOPTools_AlgoTools.OrientFacesOnShell(commonshell)
    
    shu=ShapeUpgrade_UnifySameDomain(commonshell)
    shu.Build()
    commonshell=shu.Shape()
    
    return commonshell



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
            brepgprop.SurfaceProperties(f,gpp)
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
           
            brepgprop.SurfaceProperties(f,gpp)
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



    
  
    
 
