import sys
sys.path.append('../')
from rtaa_ifc import ventilation


if __name__ == "__main__":

    filename = 'data/DCE_CDV_BAT.ifc'

    rsv=ventilation.rtaa_ventilation_study(filename)
    lgt4b=[499,705,731,757,783,809,859]
    lgt4b=[499]


    rsv.add_space_elements(lgt4b,[])#,['IfcSpace'])
    rsv.add_opening_elements([],['IfcWindow','IfcDoor'])


    rsv.set_geometries()
    rsv.run()
    rsv.display()