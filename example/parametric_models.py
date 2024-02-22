import sys
sys.path.append('../')
from rtaa_ifc.rtaa_models import *
from rtaa_ifc import solar

if __name__ == "__main__":
    lp=1.0
    hp=1.0
        
    face= ref_face(lp,hp)
    #casquette horizontale seule
    cas_h=casquette_horizontale(lp=lp,hp=hp,dhp=.1,pcg=0.1,pcd=0.1,dhm=.3,e=0.01)
    # joues rectangulaires
    gr,dr=joues_rectangulaires(lp=lp,hp=hp,ehg=0.1,ehd=0.1,dvg=.2,dvd=0.2,e=0.01)
    # joues triangulaire avec cote sup horizontal
    gt,dt=joues_triangulaires(lp=lp,hp=hp,dvg=.5,dvd=0.5)

    #casquette inclinee
    alpha=20
    cas_i=casquette_inclinee(lp=1.,hp=1.,dhp=.2,dhm=1,alpha=alpha)
    #joue triangulaire avec cote sup incline
    gti,dti=joues_triangulaires_inclinee(lp=1.,hp=1.,dhmg=1.,dhmd=1.,dhp=.2,alpha=alpha)


            
    rss=solar.rtaa_solar()
    rss.set_geometry(face,[cas_i,gti,dti])
    rss.display()
    rss.run()