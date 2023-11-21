import ifcopenshell
import sys

sys.path.append('../../IfcOpenshell/src/ifctester')

import ifctester
import ifctester.reporter
from  ifctester import ids


my_ids = ids.Ids(
    title="RTAA DOM",
    copyright="Universite de la Reunion",
    version="0.0.1",
    description="Donnees necessaire a la verification automatique par le module xxx",
    author="voivret66@gmail.com",
    date="2022-01-01",
    #purpose="Contractual requirements",
    )

location_spec = ids.Specification(
    name="Project Location",
    minOccurs=0,
    maxOccurs="unbounded",
    ifcVersion="IFC4",
    description="Localisation geographique valide du projet (lat,lon,elev)",
    instructions="Mettre à jour les coordonnées du projet (lat,lon,elev)",
    )
    
location_spec.applicability.append(ids.Entity(name="IFCSITE"))

latitude = ids.Attribute(
    name="RefLatitude",
    instructions="Le createur du projet est responsable de la bonne localisation du batiment",
    )
longitude = ids.Attribute(
    name="RefLongitude",
    instructions="Le createur du projet est responsable de la bonne localisation du batiment",
    )
altitude = ids.Attribute(
    name="RefElevation",
    instructions="Le createur du projet est responsable de la bonne localisation du batiment",
    )

location_spec.requirements.extend([latitude,longitude,altitude])

my_ids.specifications.append(location_spec)


### window
window_spec = ids.Specification(
    name="Opening porosity",
    minOccurs=0,
    maxOccurs="unbounded",
    ifcVersion="IFC4",
    description="Existence la porosite",
    instructions="test",
    )
window_spec.applicability.append(ids.Entity(name="IFCWINDOW"))

  
porosity = ids.Property(
    name='porosite',
    propertySet="Pset_rtaadom_porosite",
    datatype="IfcReal",
    instructions="sdcsct",
    )
    
window_spec.requirements.append(porosity)

my_ids.specifications.append(window_spec)


"""
door_spec = ids.Specification(
    name="Opening porosity",
    minOccurs=1,
    maxOccurs="unbounded",
    ifcVersion="IFC4",
    description="Existence la porosite",
    instructions="test",
    )
door_spec.applicability.append(ids.Entity(name="IFCDOOR"))

  
porosity = ids.Property(
    name="porosite",
    propertySet="Pset_rtaadom_porosite",
    datatype="IfcReal",
    instructions="sdcsct",
    )
door_spec.requirements.append(porosity)

my_ids.specifications.append(door_spec)
"""

if __name__ == "__main__":
    
    result = my_ids.to_xml('../tests/data/rtaa.ids')
    
    filepath =  '../tests/data/debords_casquettes_fins_triangle_pset.ifc'
    ifc=ifcopenshell.open(filepath)
    
    my_ids.validate(ifc)
    
    engine = ifctester.reporter.Console(my_ids)
    engine.report()
    #engine.to_file('./')
    #print( engine.to_string())







