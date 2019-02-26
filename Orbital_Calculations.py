
#%%
import numpy as np
import math
from numpy import cos as c
from numpy import sin as s

def orbital_cal(a, i, o, f, v1, m1, m2):  #we assume that the primary star is at z=0 with the velocity v1
    Rot=np.array([[c(f)*c(o)-c(i)*s(f)*s(o),    -c(f)*s(o)-c(i)*c(o)*s(f),  s(f)*s(i)],
                  [c(o)*s(f)+c(f)*c(i)*s(o),    c(f)*c(i)*c(o)-s(f)*s(o),   -c(f)*s(i)],
                  [s(i)*s(o),                   c(o)*s(i),                  c(i)]])
    # print('\nRot=\n',Rot)
    r0=a*np.array([[1],[0],[0]])
    # print('\nr0=\n',r0)
    r=np.matmul(Rot, r0)
    # print('\nr=\n',r)
    v0=np.array([[np.sqrt((m1**2)/(a*(m1+m2)))],[0],[0]])
    # print('\nv0= \n',v0)
    vs=(np.matmul(Rot, v0)+v1)
    # print('\nrotated velocity= \n', np.matmul(Rot,v0))
    # print('\nvs= \n', vs)
    return (r, vs)
#%%
if __name__ == "__main__":
    parameters=orbital_cal(a=1, i=math.pi/6, o=0, f=math.pi/4, v1=np.array([[0],[0],[-2]]), m1=1, m2=0.5)
    print('The position of the secondary star: \n',parameters[0], '\n\n velocity of the secondary star: \n', parameters[1])
    pass
