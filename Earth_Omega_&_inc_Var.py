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

'''here goes the omega and inclination counts we need'''
omega_count= 100            #number of omegas to try from pi/(omega_count/2) to 2*pi
inc_count= 100              #number of inclinations -1 to try from -pi/2 to pi/2 (we have 30 steps)


'''here we make lists of omega,x,vz,inc and m'''
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
x_list = get_list(3.9, 0.1, 3.9)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)
inc_list = get_list(r=-(math.pi/2), s=math.pi/(inc_count), lim=math.pi/2, rel=operator.le)
m_list = get_list(r=1, s=0.5, lim=1)
print(len(inc_list))
print(len(omega_list))
Percentage_array=[]

res = {}
with open('/Users/Behzadarbab/Exoplanet_Simulations/Earth_Omega&inc_Var/Data_x({min}-{max})_{total_count}k.txt'.format(min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as f:
    for m in m_list:
        print("m= ", m)

        for x in x_list:

            print("x= ", x)
            plt.figure(x)

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
                        sim.add(m=m, x=x, z=50, vz=v)                          #add the passing star

                        sim.integrate(round(2 * 50 / (-v), 0))                           #integration
                        a_e = sim.calculate_orbits()[0]                         #calculating orbital elements after the integration
                        res[(m, x, v, inc, o)] = (a_e.a, a_e.e)

                        '''determine the color of HZ orbits in omega-inc plot'''
                        if (a_e.a*(1+a_e.e))<1.2 and (a_e.a*(1-a_e.e))>0.8:
                            Color = 'g'
                            counter += 1
                        else:
                            Color = 'y'


                        plt.plot(o,inc,color=Color,marker='.',ms=2)

                        f.write(
                            '{m}\t{x}\t{v}\t{inc:.3f}\t{o:.3f}\t{a:.4f}\t{e:.4f}\n'.format(m=m, x=x, v=v, inc=inc, o=o,
                                                                                            a=a_e.a,
                                                                                            e=a_e.e))
                Percentage=(round(counter/(omega_count*(inc_count))*100, 2))
                # Percentage_array.append(Percentage)
                plt.xlabel('Omega (initial phase)')
                plt.ylabel('inclination')
                plt.title('X={x}, m={m}, vz={vz}({total_count}k)({P} percent)'.format(x=round(x,1),m=m, vz=v, total_count=(omega_count*(inc_count)/1000), P=Percentage))
                # plt.axes().set_aspect('equal', 'datalim')
                plt.savefig('/Users/Behzadarbab/Exoplanet_Simulations/Earth_Omega&inc_Var/X{x}m{m}v{v}({total_count}k)({P} percent).pdf'.format(x=round(x,1),m=m, v=v, total_count=(omega_count*(inc_count)/1000),
                                                                                                                                    P=Percentage), bbox_inches='tight')
