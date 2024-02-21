import sys


import numpy as np  
from collections import defaultdict
import pandas as pd
import multiprocessing


from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepBuilderAPI import (
                                    BRepBuilderAPI_Copy,
                                    BRepBuilderAPI_Transform,
                                    BRepBuilderAPI_MakeFace,
                                    BRepBuilderAPI_MakePolygon
                                    )
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.BRepGProp import brepgprop #_SurfaceProperties,brepgprop_VolumeProperties
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_Plane
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common

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
from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color,Quantity_TOC_RGB

def ref_face(lp=1.0,hp=1.0):
    # lp : largeur  
    # hp : hauteur
    p1 = gp_Pnt(0.0,0.0,0.0)
    p2 = gp_Pnt(lp,0.0,0.0)
    p3 = gp_Pnt(lp,0.0,hp)
    p4 = gp_Pnt(0.0,0.0,hp)
    points=[p1,p2,p3,p4,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    return face

# a completer avec la taille de la face    
def persienne(dist,width,thick,length,angle,step,n):
    p1 = gp_Pnt(0.0,dist-width*.5,-thick*.5)
    p2 = gp_Pnt(0.0,dist+width*.5,-thick*.5)
    p3 = gp_Pnt(0.0,dist+width*.5, thick*.5)
    p4 = gp_Pnt(0.0,dist-width*.5, thick*.5)
    points=[p1,p2,p3,p4,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    prism=BRepPrimAPI_MakePrism(face,gp_Vec(1.0,0.0,0.0),False,False).Shape()
    
    rotAxis=gp_Ax1(gp_Pnt(0.,dist,0.),gp_Dir(1.,0.,0.))
    rot= gp_Trsf()
    rot.SetRotation(rotAxis,angle)
    prism=BRepBuilderAPI_Transform(prism,rot).Shape()
    
    lshape=[prism]
    translation=gp_Trsf()
    for i in range(1,n):
        translation.SetTranslation(gp_Vec(0.0,0.0,1.0)*i*step)
        s=BRepBuilderAPI_Transform(prism,translation,True).Shape()
        lshape.append(s)
    
    return lshape
    
def casquette_horizontale(lp=1.,hp=1.,dhp=.1,pcg=0.0,pcd=0.0,dhm=.1,e=0.01):
    # create a rectangle and extrude
    Xmin= -pcd
    Xmax= lp + pcg
    Zmin= hp + dhp
    Zmax= hp + dhp + e
    p1 = gp_Pnt(Xmin,0.0,Zmin)
    p2 = gp_Pnt(Xmax,0.0,Zmin)
    p3 = gp_Pnt(Xmax,0.0,Zmax)
    p4 = gp_Pnt(Xmin,0.0,Zmax)
    points=[p1,p2,p3,p4,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    prism=BRepPrimAPI_MakePrism(face,gp_Vec(0.0,dhm,0.0),False,False).Shape()
    return prism
    

    
def joues_rectangulaires(lp=1.,hp=1.,dpg=.0,dpd=0.0,ehg=0.0,ehd=0.0,dvg=.1,dvd=0.1,e=0.01):
    # create a rectangle on XZ and extrude along Y
    Xmin= -dpd -e
    Xmax= -dpd 
    Zmin= 0.0
    Zmax= hp + ehd 
    p1 = gp_Pnt(Xmin,0.0,Zmin)
    p2 = gp_Pnt(Xmax,0.0,Zmin)
    p3 = gp_Pnt(Xmax,0.0,Zmax)
    p4 = gp_Pnt(Xmin,0.0,Zmax)
    points=[p1,p2,p3,p4,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    d=BRepPrimAPI_MakePrism(face,gp_Vec(0.0,dvd,0.0),False,False).Shape()


    Xmin= lp+ dpg 
    Xmax= lp+ dpg + e 
    Zmin= 0.0
    Zmax= hp + ehg 
    p1 = gp_Pnt(Xmin,0.0,Zmin)
    p2 = gp_Pnt(Xmax,0.0,Zmin)
    p3 = gp_Pnt(Xmax,0.0,Zmax)
    p4 = gp_Pnt(Xmin,0.0,Zmax)
    points=[p1,p2,p3,p4,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    g=BRepPrimAPI_MakePrism(face,gp_Vec(0.0,dvg,0.0),False,False).Shape()
        
    
    return d,g
    
    
def joues_triangulaires(lp=1.,hp=1.,dpg=.0,dpd=0.0,ehg=0.0,ehd=0.0,dvg=.1,dvd=0.1,e=0.01):
    #define a triangle in YZ and extrude along X
     
    p1 = gp_Pnt(-dpd,0.0,0.0)
    p2 = gp_Pnt(-dpd,0.0,hp+ehd)
    p3 = gp_Pnt(-dpd,dvd,hp+ehd)
        
    points=[p1,p2,p3,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    d=BRepPrimAPI_MakePrism(face,gp_Vec(-e,0,0.0),False,False).Shape()

    # joue gauche
    p1 = gp_Pnt(lp+dpg,0.0,0.0)
    p2 = gp_Pnt(lp+dpg,0.0,hp+ehg)
    p3 = gp_Pnt(lp+dpg,dvg,hp+ehg)
        
    points=[p1,p2,p3,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    g=BRepPrimAPI_MakePrism(face,gp_Vec(e,0,0.0),False,False).Shape()
    
    return d,g

def casquette_inclinee(lp=1.,hp=1.,dhp=0.0,pcg=0.0,pcd=0.0,dhm=.1,alpha=0.0,e=0.01):
    # create a rectangle and extrude
    Xmin= -pcd
    Xmax= lp + pcg
    Zmin= hp + dhp
    Zmax= hp + dhp + e
    p1 = gp_Pnt(Xmin,0.0,Zmin)
    p2 = gp_Pnt(Xmax,0.0,Zmin)
    p3 = gp_Pnt(Xmax,0.0,Zmax)
    p4 = gp_Pnt(Xmin,0.0,Zmax)
    points=[p1,p2,p3,p4,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    ext_vec=gp_Vec(0.0,dhm,-dhm*np.tan(np.deg2rad(alpha)))
    prism=BRepPrimAPI_MakePrism(face,ext_vec,False,False).Shape()
    return prism        
    
def joues_triangulaires_inclinee(lp=1.,hp=1.,dpg=0.0,dpd=0.0,ehg=0.0,ehd=0.0,dhmd=0.1,dhmg=0.1,dhp=0.0,alpha=45,e=0.01):
    #define a triangle in YZ and extrude along X
        
    p1 = gp_Pnt(-dpd,0.0,0.0)
    p2 = gp_Pnt(-dpd,0.0,hp+ehd+dhp)
    p3 = gp_Pnt(-dpd,dhmd,hp+ehd+dhp-dhmd*np.tan(np.deg2rad(alpha)))
        
    points=[p1,p2,p3,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    d=BRepPrimAPI_MakePrism(face,gp_Vec(-e,0,0.0),False,False).Shape()

    # joue gauche
    p1 = gp_Pnt(lp+dpg,0.0,0.0)
    p2 = gp_Pnt(lp+dpg,0.0,hp+ehg+dhp)
    p3 = gp_Pnt(lp+dpg,dhmg,hp+ehg+dhp-dhmg*np.tan(np.deg2rad(alpha)))
        
    points=[p1,p2,p3,p1]
    points.reverse()
    poly_builder = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly_builder.Add(p)
        
    poly_builder.Build()
    wire=poly_builder.Wire()
    
    face=BRepBuilderAPI_MakeFace(wire).Shape()
    g=BRepPrimAPI_MakePrism(face,gp_Vec(e,0,0.0),False,False).Shape()
    
    return d,g

