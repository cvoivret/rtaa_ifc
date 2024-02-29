import sys
sys.path.append('../')
from rtaa_ifc import solar


if __name__ == "__main__":
    filename = 'data/debords_casquettes_fins.ifc'

    rss=solar.rtaa_solar_ifc(filename)
    rss.add_building_elements([],['IfcWall','IfcSlab'])
    rss.add_solar_elements([],['IfcWindow'])
    rss.set_geometries()
    rss.run()

