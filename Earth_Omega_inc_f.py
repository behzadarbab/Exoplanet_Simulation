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
inc_count= 48              #number of inclinations to try from -pi/2 to pi/2 (we have 30 steps)
f_count= 48                #number of true anomalies from the step size to 2*pi


''''here we make lists of mass[m], Impact parameter [x], secondary star velocity [vz], 
   orbital inclination [inc], longitude of node [omega] and True anomaly [f]         '''

m_list = get_list(r=0.5, s=0.5, lim=0.5)
x_list = get_list(0.4, 1, 3)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)

inc_list = get_list(r=-((math.pi/2)-math.pi/(inc_count)), s=math.pi/(inc_count), lim=math.pi/2, rel=operator.le)
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
f_list = get_list(r=2*math.pi/f_count, s=2*math.pi/f_count, lim=2*math.pi)

'''Printing the lists lengths to see if there is any problem'''
print(inc_list)
print('length of inclination list= {inclist}'.format(inclist=len(inc_list)))
print('length of Omega list= {omegalist}'.format(omegalist=len(omega_list)))
print('length of True anomaly list= {anolist}'.format(anolist=len(f_list)))


Percentage_array=[]

'''calculating the dimensions of the figure (subplots in height and width)'''
subnumber = math.ceil(math.sqrt(len(inc_list)))
print('number of rows and columns of subplots: {sub}'.format(sub=subnumber))

res = {}
with open('/Users/Behzadarbab/Exoplanet_Simulations/Earth_m_Omega_inc_f/chances_m({mmin}-{mmax})_x({min}-{max})_{total_count}k.txt'.format(mmin=min(m_list), mmax=max(m_list), min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as g:
    g.write('x\t%\n')
    with open('/Users/Behzadarbab/Exoplanet_Simulations/Earth_m_Omega_inc_f/Data_m({mmin}-{mmax})_x({min}-{max})_{total_count}k.txt'.format(mmin=min(m_list), mmax=max(m_list), min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as h:
       for m in m_list:
            print("m= ", m)

            for x in x_list:

                print("x= ", x)

                plt.figure(figsize=(7*subnumber,7*subnumber))

                counter = 0  # this is the counter for calculating percentage

                for v in vz_list:

                    print("vz= ", v)


                    for inc in inc_list:

                        plotnumber= int(round((inc+math.pi/2)*inc_count/math.pi, 0))
                        print('plot number: {pl}'.format(pl=plotnumber))

                        plt.subplot(subnumber,subnumber, plotnumber)

                        inc_name= (inc+math.pi/2)*inc_count/(math.pi)
                        print("\t inc: {Q:.0f}".format(Q=inc_name))



                        for o in omega_list:

                            for f in f_list:

                                # print('\t\t Omega: {O:.0f}'.format(O=o*omega_count/(2*math.pi)))

                                sim = rebound.Simulation()                              #the simulation starts here

                                sim.add(m=1)                                            #add Sun
                                sim.add(m=3e-6, a=1, Omega=o, inc=inc, f=f)             #add Earth
                                sim.add(m=m, x=x, z=50, vz=v)                           #add the passing star

                                sim.integrate(round(2 * 50 / (-v), 0))                  #integration

                                a_e = sim.calculate_orbits()[0]                         #calculating orbital elements after the integration
                                res[(m, x, v, inc, o, f)] = (a_e.a, a_e.e)

                                '''determine the color of HZ orbits in omega-inc plot'''
                                if (a_e.a*(1+a_e.e))<1.2 and (a_e.a*(1-a_e.e))>0.8:
                                    Color = '0'
                                    counter += 1
                                else:
                                    Color = '0.75'


                                plt.plot(o, f, color=Color, marker='s', ms=6)


                                h.write(
                                    '{m}\t{x}\t{v}\t{inc:.3f}\t{o:.3f}\t{f:.3f}\t{a:.3f}\t{e:.3f}\n'.format(m=m, x=x, v=v, inc=inc, o=o, f=f,
                                                                                            a=a_e.a,
                                                                                            e=a_e.e))

                        plt.xlabel('Omega (initial phase)')
                        plt.ylabel('True Initial Anomaly')
                        plt.title('inc={inc}[rad]'.format(inc=round(inc, 3)))
                    plt.tight_layout()
                    plt.subplots_adjust(top=0.95)
                    Percentage = (round((counter / (inc_count * omega_count * (f_count))) * 100, 2))
                    plt.suptitle('X={x}, m={m}, vz={vz},({total_count}k) Percentage={P}'.format(x=round(x, 2), m=m, vz=v,total_count=(omega_count*inc_count*f_count), P=Percentage), fontsize=25)
                    plt.savefig('/Users/Behzadarbab/Exoplanet_Simulations/Earth_Omega_inc_f/X{x}m{m}v{v}({total_count}k).pdf'.format(x=round(x,2),m=m, v=v, total_count=(omega_count*(inc_count)*f_count/1000)
                                                                                                                                        ), bbox_inches='tight')
                    g.write('{x}\t{Percentage}\n'.format(x=x, Percentage=Percentage))
                plt.close()