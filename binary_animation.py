import rebound
import math
import operator
import matplotlib.pyplot as plt
import numpy as np
from numpy import cos as c
from numpy import sin as s

'''we define listing here'''
def get_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(r)
        r += s
    return my_list

def orbital_cal(a, i, o, f, vz, m1, m2):  #we assume that the primary star is at z=0 with the velocity -vz
    Rot=np.array([[c(f)*c(o)-c(i)*s(f)*s(o),    -c(f)*s(o)-c(i)*c(o)*s(f),  s(f)*s(i)],
                  [c(o)*s(f)+c(f)*c(i)*s(o),    c(f)*c(i)*c(o)-s(f)*s(o),   -c(f)*s(i)],
                  [s(i)*s(o),                   c(o)*s(i),                  c(i)]])
    print('Rot=',Rot)
    r0=a*np.array([[1],[0],[0]])
    print('r0=',r0)
    r=np.matmul(Rot, r0)
    print('r=',r)
    v0=np.array([[np.sqrt((m1**2)/(a*(m1+m2)))],[0],[0]])
    print('v0=',v0)
    vs=np.matmul(Rot, v0)-np.array([[0],[0],[-vz]])
    print('vs=', vs)
    return (r, vs)

'''here goes the Omega, inclination and true anomaly counts we need'''
omega_count= 4            #number of omegas to try from pi/(omega_count/2) to pi (we dont use pi to 2pi because it is symmetric
inc_count= 4               #number of inclinations to try from -pi/2 to pi/2 (we have 30 steps)
f_count= 4                #number of true anomalies from the step size to 2*pi


'''here we make lists of omega,x,vz,inc and m'''
omega_list = get_list(-math.pi, 2* math.pi/omega_count, math.pi-math.pi/omega_count, rel=operator.le)
x_list = get_list(5, 1, 5)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)
inc_list = get_list(r=0, s=math.pi/(inc_count), lim=((math.pi+0.01/2)), rel=operator.le)
m_list = get_list(r=1, s=0.5, lim=1)
f_list = get_list(r=-math.pi, s=2*math.pi/f_count, lim=math.pi-math.pi/f_count)

'''Printing the lists lengths to see if there is any problem'''
print(inc_list)
print('length of inclination list= {inclist}'.format(inclist=len(inc_list)))
print('length of Omega list= {omegalist}'.format(omegalist=len(omega_list)))
print('length of True anomaly list= {anolist}'.format(anolist=len(f_list)))

# Percentage_array=[]

'''calculating the dimensions of the figure (subplots in height and width)'''

res = {}

for m in m_list:
    print("m= ", m)

    for x in x_list:

        print("x= ", x)

        for v in vz_list:
            print("vz= ", v)
            counter = 0
            for inc in inc_list:
                inc_name= inc*inc_count/(math.pi)
                print("\t inc: {Q:.0f}".format(Q=inc_name))

                for o in omega_list:
                    for f in f_list:
                        print ('f= {f}'.format(f=f))
                        sim = rebound.Simulation()  # the simulation starts here
                        sim.move_to_com()
                        sim.add(m=1)  # add Sun
                        sim.add(m=0, a=1, Omega=o, inc=inc, f=f)  # add Earth

                        (r,vs)=orbital_cal(a=1, i=(np.pi)/4, o=0, f=0, vz=1, m1=1, m2=0.5)
                        v1= (0.5/1)*(np.array([[0],[0],[-100]])-vs)
                        sim.add(m=1, x=x, y=0, z=50, vy=-math.sqrt(1/20),vz=v1)  # add the passing star


                        sim.add(primary=sim.particles[2], m=1, a=3)


                        for tt in range(100):
                            print('t= {t}'.format(t=tt))

                            sim.integrate(round((tt)/((-v)), 0))                           #integration

                            fig, ax = plt.subplots(figsize=(8, 8))
                            ps = sim.particles
                            sim.move_to_com()

                            # manually set plot boundaries
                            lim = 20
                            ax.set_xlim([-lim, lim])
                            ax.set_ylim([-lim, lim])

                            # plot the stars and planets with separate symbols
                            linewidth = 1.


                            ax.scatter(ps[0].x, ps[0].y, s=40 * linewidth, marker='*', facecolor='black', zorder=3)
                            ax.scatter(ps[2].x, ps[2].y, s=40*(ps[2].m)/(ps[0].m) * linewidth, marker='*', facecolor='black', zorder=3)
                            ax.scatter(ps[3].x, ps[3].y, s=40*(ps[3].m)/(ps[0].m) * linewidth, marker='*', facecolor='black', zorder=3)

                            ax.scatter(ps[1].x, ps[1].y, s=10 * linewidth, facecolor='blue', zorder=3)

                            # Now individually plot orbit trails with appropriate orbit

                            from rebound.plotting import fading_line

                            planet = ps[1]  # circumbinary planet, use default jacobi coordinates
                            o1 = np.array(planet.sample_orbit())
                            lc = fading_line(o1[:, 0], o1[:, 1], linewidth=linewidth)
                            ax.add_collection(lc)

                            secondary = ps[3]  # planet in orbit around B, assign it as primary
                            o1 = np.array(secondary.sample_orbit(primary=ps[2]))
                            lc = fading_line(o1[:, 0], o1[:, 1], linewidth=linewidth)
                            ax.add_collection(lc)
                            # fig = rebound.OrbitPlot(sim, trails=True, slices=True, color=True, periastron=False,
                            #                         unitlabel="[AU]", lim=30., limz=30)

                            plt.suptitle('X={x} m={m} vz={vz} inc={inc} Omega={Omega} Anomaly={f}'.format(x=round(x,1), m=m, vz=v, inc=inc, Omega=o, f=f))
                            plt.savefig(
                                '/Users/atefeh-behzad/Exoplanet_Simulations/binary_encounter/tests/1X{x}m{m}v{v}inc{inc}Omega{o}Anomaly{f}t{t}.png'.format(
                                    x=round(x, 1), m=m, v=v, inc=round(inc, 2), o=round(o, 2), f=round(f, 2), t=(tt+25)), dpi=100)
                            # plt.show()
                            plt.close()