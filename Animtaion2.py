import rebound
import math
import operator
import matplotlib.pyplot as plt



'''we define listing here'''
def get_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(r)
        r += s
    return my_list

'''here goes the Omega, inclination and true anomaly counts we need'''
omega_count= 48            #number of omegas to try from pi/(omega_count/2) to pi (we dont use pi to 2pi because it is symmetric
inc_count= 1              #number of inclinations to try from -pi/2 to pi/2 (we have 30 steps)
f_count= 48                #number of true anomalies from the step size to 2*pi


'''here we make lists of omega,x,vz,inc and m'''
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
x_list = get_list(0.5, 1, 0.5)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)
inc_list = get_list(r=-(math.pi/2)+9*math.pi/(48), s=math.pi/(48), lim=-((math.pi/2)-9.02*math.pi/(48)), rel=operator.le)
m_list = get_list(r=1, s=0.5, lim=1)
f_list = get_list(r=2*math.pi/f_count, s=2*math.pi/f_count, lim=2*math.pi)

'''Printing the lists lengths to see if there is any problem'''
print(inc_list)
print('length of inclination list= {inclist}'.format(inclist=len(inc_list)))
print('length of Omega list= {omegalist}'.format(omegalist=len(omega_list)))
print('length of True anomaly list= {anolist}'.format(anolist=len(f_list)))


# Percentage_array=[]

'''calculating the dimensions of the figure (subplots in height and width)'''
# subnumber = math.ceil(math.sqrt(len(inc_list)))
# print('number of rows and columns of subplots: {sub}'.format(sub=subnumber))

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

                        sim.add(m=1)  # add Sun
                        sim.add(m=3e-6, a=1, Omega=o, inc=inc, f=f)  # add Earth
                        sim.add(m=m, x=x, z=50, vz=v)  # add the passing star

                        sim.move_to_com()
                        for tt in range(50):
                            print('t= {t}'.format(t=tt))

                            # print('\t\t Omega: {O:.0f}'.format(O=o*omega_count/(2*math.pi)))


                            sim.integrate(round(2 * 50*(tt+25) / ((-v)*100), 0))                           #integration


                            fig = rebound.OrbitPlot(sim, trails=True, slices=True, color=True,periastron=True, unitlabel="[AU]", lim=5., limz=5)
                            plt.suptitle('X={x} m={m} vz={vz} inc={inc} Omega={Omega} Anomaly={f}'.format(x=round(x,1), m=m, vz=v, inc=inc, Omega=o, f=f))
                            plt.savefig(
                                '/Users/Behzadarbab/Exoplanet_Simulations/Earth_Omega&inc_Var_Animation/6/X{x}m{m}v{v}inc{inc}Omega{o}Anomaly{f}t{t}.png'.format(
                                    x=round(x, 1), m=m, v=v, inc=round(inc, 2), o=round(o, 2), f=round(f, 2), t=(tt+25)), dpi=100)
                            plt.close()