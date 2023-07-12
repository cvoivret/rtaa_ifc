import numpy as np
import scipy.integrate as spi

#theta : 0 Pi
#phi : 0 2Pi

n=50
t= np.linspace(0,np.pi,n,endpoint=True)
p= np.linspace(0,2*np.pi,2*n,endpoint=True)

tp=np.meshgrid(t,p)

ds=lambda x,y: np.sin(x)

result, error = spi.nquad(ds, [[0, np.pi],[0, 2*np.pi]])

surf = np.empty((t.shape[0],p.shape[0]))
for i in range(t.shape[0]-1):
    for j in range(p.shape[0]-1):
        res,err = spi.nquad(ds, [[t[i], t[i+1]],[p[j], p[j+1]]])
        surf[i,j]=res
        


