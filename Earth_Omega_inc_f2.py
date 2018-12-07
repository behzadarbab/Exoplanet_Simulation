"""                                             Briefing about the simulation
This program simulates a three body system; Sun, Earth, a passing star. The passing star starts from 50 AU above the
Solar system plane. We throw the star with a speed (or we can change the speed or give a list of speeds), then after a
time that it reaches around 50 AU below the solar plane, we stop the single simulation and check if the planet orbit
reaches beyond the HZ aroound sun (0.8AU-1.2AU). Then we do this simulation for different inclinations and directions
of the planet orbit and we change the place of the planet on the orbit. After this, we can make a 3D plot (a grand plot
with many subplots in 2D) and if in a simulation, after the flyby, the planet stays in HZ, we plot it in black, otherwise
in gray."""



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
omega_count= 10            #number of omegas to try from pi/(omega_count/2) to pi (we dont use pi to 2pi because it is symmetric
inc_count= 10              #number of inclinations to try from -pi/2 to pi/2 (we have 30 steps)
f_count= 10                #number of true anomalies from the step size to 2*pi


'''here we make lists of omega,x,vz,inc and m'''
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
x_list = get_list(0.5, 0.01, 10)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)
inc_list = get_list(r=-((math.pi/2)-math.pi/(inc_count)), s=math.pi/(inc_count), lim=math.pi/2, rel=operator.le)
m_list = get_list(r=1, s=0.5, lim=1)
f_list = get_list(r=2*math.pi/f_count, s=2*math.pi/f_count, lim=2*math.pi)

'''Printing the lists lengths to see if there is any problem'''
print(inc_list)
print('length of inclination list= {inclist}'.format(inclist=len(inc_list)))
print('length of Omega list= {omegalist}'.format(omegalist=len(omega_list)))
print('length of True anomaly list= {anolist}'.format(anolist=len(f_list)))


Percentage_array=[]

'''calculating the dimensions of the figure (subplots in height and width)'''
# subnumber = math.ceil(math.sqrt(len(inc_list)))
# print('number of rows and columns of subplots: {sub}'.format(sub=subnumber))

res = {}
with open('/Users/atefeh-behzad/Exoplanet_Simulations/Earth_Omega_inc_f_2/chances_x({min:.2f}-{max:.2f})_{total_count}k.txt'.format(min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as g:
    g.write('x\t%\n')
    with open('/Users/atefeh-behzad/Exoplanet_Simulations/Earth_Omega_inc_f_2/Data_x({min:.2f}-{max:.2f})_{total_count}k.txt'.format(min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as h:
        h.write('StarM\tStarx\tStarv\tPiinc\tPiOmega\tPif\tPfinc\tPfOmega\tPff\taf\tef\n')
        for m in m_list:
            print("m= ", m)

            for x in x_list:

                print("x= ", x)

                # plt.figure(figsize=(7*8,7*6))

                counter = 0  # this is the counter for calculating percentage

                for v in vz_list:

                    print("vz= ", v)


                    for inc in inc_list:

                        plotnumber= int(round((inc+math.pi/2)*inc_count/math.pi, 0))
                        print('\nplot number: {pl}'.format(pl=plotnumber))

                        # plt.subplot(6,8, plotnumber)

                        inc_name= (inc+math.pi/2)*inc_count/(math.pi)
                        # print("\t inc: {Q:.0f}".format(Q=inc_name))



                        for o in omega_list:
                            print('\nx={x}, Inc={inc:.0f}, Omega: {o:.0f}'.format(x=x, inc=inc_name,o=omega_count*o/(2*math.pi)))
                            print("", end="                         Anomaly # ")
                            for f in f_list:
                                print("{r:.0f}".format(r=f_count*f/(2*math.pi)), end=" ")

                                # print('\t\t Omega: {O:.0f}'.format(O=o*omega_count/(2*math.pi)))

                                sim = rebound.Simulation()                              #the simulation starts here

                                sim.add(m=1)                                            #add Sun
                                sim.add(m=3e-6, a=1, Omega=o, inc=inc, f=f)             #add Earth
                                sim.add(m=m, x=x, z=50, vz=v)                           #add the passing star
                                a_c = sim.calculate_orbits()[0]
                                print('initial1: 1.00, 0.00, {inc:.2f}, {o:.2f}, {f:.2f}\n'.format(o=o, inc=inc, f=f))
                                print('  initial2: {a:.2f}, {e:.2f}, {inc:.2f}, {Omega:.2f}, {f:.2f} \n'.format(a=a_c.a, e=a_c.e, inc=a_c.inc, Omega=a_c.Omega, f=a_c.f))
                                sim.integrate(round(2 * 50 / (-v), 0))                  #integration

                                a_e = sim.calculate_orbits()[0]                         #calculating orbital elements of the planet after the integration
                                res[(m, x, v, inc, o, f)] = (a_e.a, a_e.e)
                                print('  final:    {a:.2f}, {e:.2f}, {inc:.2f}, {Omega:.2f}, {f:.2f} \n'.format(a=a_e.a,
                                                                                                               e=a_e.e,
                                                                                                               inc=a_e.inc,
                                                                                                               Omega=a_e.Omega,
                                                                                                               f=a_e.f))
                                '''determine the color of HZ orbits in omega-inc plot'''
                                if (a_e.a*(1+a_e.e))<1.2 and (a_e.a*(1-a_e.e))>0.8:
                                    Color = '0'
                                    counter += 1
                                else:
                                    Color = '0.75'


                                plt.plot(o, f, color=Color, marker='s', ms=6)


                                h.write(
                                    '{m}\t{x:.2f}\t{v}\t{inc:.3f}\t{o:.3f}\t{f:.3f}\t{finc:.3f}\t{fo:.3f}\t{ff:.3f}\t{a:.3f}\t{e:.3f}\n'.format(m=m, x=x, v=v, inc=inc, o=o, f=f, finc=a_e.inc, fo=a_e.Omega, ff=a_e.f,
                                                                                            a=a_e.a,
                                                                                            e=a_e.e))

                        plt.xlabel('Omega (initial phase)')
                        plt.ylabel('True Initial Anomaly')
                        plt.title('inc={inc}[rad]'.format(inc=round(inc, 3)))
                    plt.tight_layout()
                    plt.subplots_adjust(top=0.95)
                    Percentage = (round((counter / (inc_count * omega_count * (f_count))) * 100, 2))
                    plt.suptitle('X={x}, m={m}, vz={vz},({total_count}k) Percentage={P}'.format(x=round(x, 2), m=m, vz=v,total_count=(omega_count*inc_count*f_count), P=Percentage), fontsize=25)
                    plt.savefig('/Users/atefeh-behzad/Exoplanet_Simulations/Earth_Omega_inc_f_2/X{x}m{m}v{v}({total_count}k).pdf'.format(x=round(x,2),m=m, v=v, total_count=(omega_count*(inc_count)*f_count/1000)
                                                                                                                                        ), bbox_inches='tight')
                    g.write('{x}\t{Percentage}\n'.format(x=x, Percentage=Percentage))
                plt.close()