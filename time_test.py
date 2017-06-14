import rebound
import numpy as np
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

'''here goes the omega and inclination counts we need'''
omega_count= 50            #number of omegas to try from pi/(omega_count/2) to 2*pi
inc_count= 50              #number of inclinations -1 to try from -pi/2 to pi/2 (we have 30 steps)


'''here we make lists of omega,x,vz,inc and m'''
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
x_list = get_list(5, 1, 5)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)
inc_list = get_list(r=-(math.pi/2), s=math.pi/(inc_count), lim=math.pi/2, rel=operator.le)
m_list = get_list(r=1, s=0.5, lim=1)
print(len(inc_list))
print(len(omega_list))
Percentage_array=[]
tmin, tmax = 1, 5
res = {}
with open('/Users/Behzadarbab/Exoplanet_Simulations/Earth_Omega&inc_Var/Data_x({min}-{max})_{total_count}k_timetest.txt'.format(min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                              total_count=omega_count*(inc_count)/1000), 'w+') as f:
    plt.figure(1)

    for t in range(1,26):
        print("t= ", t)
        plt.subplot(5,5,t)
        for m in m_list:
            print("m= ", m)

            for x in x_list:

                print("x= ", x)
                # plt.figure(x)

                for v in vz_list:
                    print("vz= ", v)
                    counter = 0
                    for inc in inc_list:
                        inc_name= inc*inc_count/(math.pi)
                        print("\t inc: {Q:.0f}".format(Q=inc_name))

                        for o in omega_list:

                            # print('\t\t Omega: {O:.0f}'.format(O=o*omega_count/(2*math.pi)))

                            sim = rebound.Simulation()                              #the simulation starts here

                            sim.add(m=1)                                            #add Sun
                            sim.add(m=3e-6, a=1, Omega=o, inc=inc)                  #add Earth
                            sim.add(m=m, x=x, z=40*t, vz=v)                          #add the passing star

                            sim.integrate(round(2 * 40*t / (-v), 0))                           #integration
                            a_e = sim.calculate_orbits()[0]                         #calculating orbital elements after the integration
                            res[(m, x, v, inc, o)] = (a_e.a, a_e.e)

                            '''determine the color of HZ orbits in omega-inc plot'''
                            if (a_e.a*(1+a_e.e))<1.2 and (a_e.a*(1-a_e.e))>0.8:
                                Color = 'g'
                                counter += 1
                            else:
                                Color = 'y'


                            plt.plot(o,inc,color=Color,marker=',',ms=1)

                            f.write(
                                '{m}\t{x}\t{v}\t{inc:.3f}\t{o:.3f}\t{a:.4f}\t{e:.4f}\n'.format(m=m, x=x, v=v, inc=inc, o=o,
                                                                                            a=a_e.a,
                                                                                            e=a_e.e))
                    Percentage=(round(counter/(omega_count*(inc_count))*100, 2))
                    # Percentage_array.append(Percentage)
                # plt.xlabel('Omega (initial phase)', fontsize='3', y=2)
                # plt.ylabel('inclination', fontsize='3', x=2)
                plt.xticks(fontsize='3')
                plt.yticks(fontsize='3')
                plt.title('t={t}, X={x}, ({total_count}k)({P} percent)'.format(t=t, x=round(x,1), total_count=(omega_count*(inc_count)/1000), P=round(Percentage,2)), fontsize='3', y=0.98)
                # plt.axes().set_aspect('equal', 'datalim')



plt.tight_layout()
plt.savefig('/Users/Behzadarbab/Exoplanet_Simulations/Earth_Omega&inc_Var/X{x}m{m}v{v}({total_count}k)_time({tmin}-{tmax}).pdf'.format(x=round(x,1),tmin=tmin, tmax=tmax, m=m, v=v, total_count=(omega_count*(inc_count)/1000),
                                                                                                                                    ))
