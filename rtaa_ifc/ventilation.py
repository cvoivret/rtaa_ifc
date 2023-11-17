from collections import defaultdict
from itertools import compress,tee,combinations,product
import os

#from collections import Counter

import ifcopenshell
from ifcopenshell.geom import create_shape
import pandas as pd
from  project import project_location
from geom import   (
                    ifcelement_as_solid,
                    fuse_listOfShape,
                    shapes_as_solids,
                    get_external_shell,
                    exterior_wall_normal_and_plane,
                    biggestface_along_vector,
                    setting
                    )
from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color,Quantity_TOC_RGB

from OCC.Core.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from OCC.Core.gp import gp_Pnt,gp_Dir,gp_Vec,gp_Pln,gp_Lin,gp_Trsf,gp_Ax3,gp_Ax1
from OCC.Core.gp import gp_Pnt2d,gp_Dir2d,gp_Vec2d,gp_Lin2d

from OCC.Core.Geom import Geom_Plane,Geom_Line

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut,BRepAlgoAPI_Fuse
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties,brepgprop_VolumeProperties
from OCC.Core.GProp import GProp_GProps

from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.ProjLib import projlib_Project
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf

from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve
from OCC.Core.Geom2d import Geom2d_Line
from OCC.Core.GeomAPI import geomapi_To3d
from OCC.Core.ElCLib import elclib_Parameter,elclib_To3d,elclib_Value


from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer,TopoDS_Iterator

#from OCC.Core.GeomLProp import GeomLProp_SLProps

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)
          
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

class space_ventilation_data:
    def __init__(self,ifcspace):
        self._space=ifcspace
        self._win_by_wall=defaultdict(list)
        self._wall_faces =defaultdict(list)
        self._win_faces  =defaultdict(list)
        

    
    def extract_faces(self,ifcspace,window_by_wall,wall_shapes,windows_shapes):
        """Build faces that reprensents the opening. Make the link 
        between spaces and wall. Construct part of wall that belong
        to a space and also the window that belong to a space. By 
        geometry only (expected to work even if the ifc file is minimal).

        Parameters
        ----------
        ifcspace : IfcSpace
            The space are interested in
        window_by_wall : dict( int : int) 
            Associate a list of window id for each wall id in the model
        wall_shapes : dict(int:Shape)
            Associate a shape to each wall id
        window_shapes : dict(int:Shape)
            Associate a shape to each window id
        
        Returns
        -------
        None
            
        """ 
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
        for  i,j in combinations(range(0,len(lax3)),r=2):
            coplanar= lax3[i].IsCoplanar(lax3[j],1.e-5,1.e-5)
            if coplanar:
                lcoplanar[i].append(j)
                #print(i,' ',j,' ')#,lcoplanar[i])
        newfaces=[]
        toremove=[]
        for f,indices in zip(faces,lcoplanar):
            if len(indices)>0:
                los=[f]
                [los.append(lface[i]) for i in indices]
                toremove.extend(indices)
                #print(los)
                new=fuse_listOfShape(los)
                #new=list(TopologyExplorer(new).faces())
                #print('new ',new)
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
                            #print(wall_id,' ',self._wall_faces[wall_id])
                            self._win_faces[win_id].extend(bigfaces)
                            #svd.update(wall_id,win_id,f,bigfaces)
    
    
    def opening_ratio(self,projloc):
        """Compute the opening ratio indicator defined in RTAA. 

        Parameters
        ----------
        projloc : project_location
            To know the face orientations
        
        Returns
        -------
        opening_data : pd.DataFrame   
            Data about spaces and associated openings
            
        """  
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
        
        data={'nom_espace':space_name,
              'id_espace':space_id,
              'id_ouverture':win_ids,
              'rtaaporo_ouverture':poros,
              'aire_ouverture':areas}
        win_data=pd.DataFrame(data)
        self._win_data=win_data
        
        
        
        wall_ids=list(self._win_by_wall.keys())
        space_name= [self._space.Name]*len(wall_ids)
        space_id= [self._space.id()]*len(wall_ids)
        wall_areas =[ wall_area_total[wall_id] for wall_id in  wall_ids]
        win_areas =[ win_area_by_wall[wall_id] for wall_id in  wall_ids]
        wall_angle=[wall_angle[wall_id] for wall_id in wall_ids]
        wall_orient=[ wall_orient[wall_id] for wall_id in wall_ids]
        data={'nom_espace':space_name,
              'id_espace':space_id,
              'id_paroi':wall_ids,
              'aire_paroi':wall_areas,
              'aire_ouverture':win_areas,
              'angleNord_paroi':wall_angle,
              'secteurOrient_paroi':wall_orient}
        wall_data=pd.DataFrame(data)
        wall_data.sort_values(by=['aire_paroi'],ascending=False,inplace=True)
        self._wall_data=wall_data
        
        self._opening_ratio=wall_data['aire_ouverture'].sum()/wall_data['aire_paroi'][0]
        
        
    def sweeping(self,windows_shapes):
        """Compute the sweeping indicator defined in RTAA. Distance between 
        two opening in the same room

        Parameters
        ----------
        windows_shapes : dict( Shape)
            Shapes of all windows (some kind of caching )
        
        Returns
        -------
        sweep_data : pd.DataFrame   
            Data about windows distance in the same space
            
        """         # compute and analyze intersection face and curve linking two windows
        
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
                length_between_mc[(win_id1,win_id2)]=length
        
        
        temp=list(length_between_mc.keys())
        win_id1 = [t[0] for t in temp]
        win_id2 = [t[1] for t in temp]
        l= list(length_between_mc.values())
        space_name= [self._space.Name]*len(temp)
        space_id= [self._space.id()]*len(temp)
        halfdiag = [half_diag]*len(temp)
        data={'nom_espace':space_name,
              'id_espace':space_id,
              'id_win1':win_id1,
              'id_win2':win_id2,
              'distance':l,
              'demidiag':halfdiag,
              }
        sweep_data=pd.DataFrame(data)
        self._sweep_data = sweep_data
        
        print(sweep_data)
        """
        print(' Sweeping (must be >1) : ' )
        for k,v in length_between_mc.items():
            print('     ',k[0],' ',k[1],' ',v)
        print(length_between_mc)
        """
        #self._sweeping_data=length_between_mc
        
        return sweep_data
    
      
    def info(self):
        print('Space Name ',self._space.Name)
        print('Space Id   ',self._space.id())
        print('*** Windows id by wall id of this space')
        for wall_id,win_list in self._win_by_wall.items():
            print("     wall id ",wall_id)
            for w in win_list:
                print("         win id ",w)
       
        print('\n')  


    


class rtaa_ventilation_study:
    def __init__(self,ifcfilename):
        setting=ifcopenshell.geom.settings()
        setting.set(setting.USE_PYTHON_OPENCASCADE, True)
        
        self._workingdir = os.path.dirname(ifcfilename)
        
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
        print(" Number of space for analysis: ",len(self._space_elements.keys()))
    
    def remove_spaces_elements(self,ids=[]):
        self._remove_elements(ids,types,'spaces')
            
        
    def add_opening_elements(self,ids=[],types=[]):
        self._add_elements(ids,types,'openings')
        print(" Number of elements in building: ",len(self._opening_elements.keys()))
    
    def remove_openings_elements(self,ids=[]):
        self._remove_elements(ids,types,'openings')    
                
    def remove_overlapping_spaces(self):
        """Not finished funtion, some space can have the same extrusion bases but 
        with different heigth

        Parameters
        ----------
        None
        
        Returns
        -------
        None
            
        """         #print(self._space_elements)
        sorted_el = dict(sorted(self._space_elements.items()))
        #print(sorted_el)
        shapes = [ ifcelement_as_solid(el) for el in sorted_el.values()]
        combishapes=itertools.combinations(shapes,2)
        combiid = itertools.combinations(sorted_el.keys(),2)
        map=defaultdict(list)
        for (s1,s2),(id1,id2) in zip(combishapes,combiid):
            s=BRepAlgoAPI_Common(s1,s2).Shape()
            
            solids=list(TopologyExplorer(s).solids())
            if( len(solids)>0):
                # identical spaces
                cut12=BRepAlgoAPI_Cut(s1,s2).Shape()
                cut12list=list(TopologyExplorer(cut12).solids())
                
                cut21=BRepAlgoAPI_Cut(s2,s1).Shape()
                cut21list=list(TopologyExplorer(cut21).solids())
                
                cutunion = BRepAlgoAPI_Fuse(cut12,cut21).Shape()
                cutcut=BRepAlgoAPI_Cut(s,cutunion).Shape()
                cutcutlist=list(TopologyExplorer(cutcut).solids())
                print(cutcutlist)
                print( len(cut12list),' ',len(cut21list),' ', len(cutcutlist))
                
                if( (len(cut12list)==0 )& (len(cut21list)==0)):
                    #identical space
                    map[id1].extend()
                elif ( (len(cut12list)==1 )& (len(cut21list)==0)):
                    # 1 include 2
                    # could test if the inclusion is on the same profile
                    # just a question of extrusion height
                    map[id1].append(id2)
                elif (( len(cut12list)==0 )& (len(cut21list)==1)):
                    # 2 include 1
                    map[id2].append(id1)
                
        #print(map)            
        G=nx.Graph()
        for k,v in map.items():
            # compute union of spaces
            #compare to
            print(k,' ',v)
            for vertex in v:
                G.add_edge(k,vertex)
            
        plt.show()
        
        
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
        """Display wall, spaces and windows faces used in computation

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
        """Generate a spreadsheet with data and results, both for analysis
           verification of results

        Parameters
        ----------
        None
        
        Returns
        -------
        None
            
        """        
        filename=os.path.join(self._workingdir,'ventilation_results.xlsx')
        with pd.ExcelWriter(filename,mode='w') as writer:  
            
            lwin=[r._win_data for r in self._results if hasattr(r,'_win_data')]
            win=pd.DataFrame()
            if(len(lwin)>0):
                win  = pd.concat(lwin)
                win_name=[ self._ifc_file.by_id(w_id).Name for w_id in win['id_ouverture']]
                #print(win_name)
                #print(win)
                win.insert(2,'nom_fenetre',win_name)
                #print(win)
                win.reset_index(drop=True,inplace=True)
                print(win.head())
                win.to_excel(writer, sheet_name='fenetres')
            
            lwall=[r._wall_data for r in self._results if hasattr(r,'_wall_data')]
            wall=pd.DataFrame()
            if(len(lwall)>0):
                wall = pd.concat(lwall)
                wall.reset_index(drop=True,inplace=True)
                wall.to_excel(writer, sheet_name='parois')
                
            lsweep=[r._sweep_data for r in self._results if hasattr(r,'_sweep_data')]
            sweep=pd.DataFrame()
            if len(lsweep)>0:
                sweep= pd.concat(lsweep)
                sweep.reset_index(drop=True,inplace=True)
                sweep.to_excel(writer, sheet_name='balayage')

        return win,wall,sweep
           
    
    def run(self):
        self._results=[]
        
        #self.remove_overlapping_spaces()
        #sdcds
        
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
            
        win,wall,sweep=self.export_raw()
        



if __name__ == "__main__":
    
    filename = '../tests/data/DCE_CDV_BAT.ifc'

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
