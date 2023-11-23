import sys
sys.path.append('../')
from rtaa_ifc import solar



filename = 'data/debords_casquettes_fins.ifc'

rss=solar.rtaa_solar_study(filename)
rss.add_building_elements([],['IfcWall','IfcSlab'])
rss.add_solar_elements([],['IfcWindow'])
rss.set_geometries()
rss.run()

