import ifcopenshell
import sys

sys.path.append('C:/Users/cvoivret/source/IfcOpenshell/src/ifctester')

import ifctester
import ifctester.reporter
from  ifctester import ids

def rtaa_solar_ids():

    rtaa_ids = ids.Ids(
        title="RTAA DOM",
        copyright="Universite de la Reunion",
        version="0.0.1",
        description="Donnees necessaire a la verification automatique par le module xxx",
        author="voivret66@gmail.com",
        date="2022-01-01",
        #purpose="Contractual requirements",
        )
    ### project location
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

    rtaa_ids.specifications.append(location_spec)

    return rtaa_ids
    
def rtaa_ventilation_ids():
    
    rtaa_ids = rtaa_solar_ids()

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
        name="porosite",
        propertySet="Pset_rtaadom_porosite",
        datatype="IfcReal",
        #instructions="sdcsct",
        )
        
    window_spec.requirements.append(porosity)

    rtaa_ids.specifications.append(window_spec)


    ### door
    door_spec = ids.Specification(
        name="Opening porosity",
        minOccurs=0,
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

    rtaa_ids.specifications.append(door_spec)
    
    return rtaa_ids


def check_solar_ids(ifc):
    rtaa_ids=rtaa_solar_ids()
    rtaa_ids.validate(ifc)
    engine = ifctester.reporter.Console(rtaa_ids)
    engine.report()
    
    return 
    
def check_ventilation_ids(ifc):
    rtaa_ids=rtaa_ventilation_ids()
    rtaa_ids.validate(ifc)
    engine = ifctester.reporter.Console(rtaa_ids)
    engine.report()
    
    return 





