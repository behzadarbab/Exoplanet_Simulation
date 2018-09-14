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
omega_count= 12            #number of omegas to try from pi/(omega_count/2) to pi (we dont use pi to 2pi because it is symmetric
inc_count= 12              #number of inclinations to try from -pi/2 to pi/2 (we have 30 steps)
f_count= 12                #number of true anomalies from the step size to 2*pi


'''here we make lists of mass[m], Impact parameter [x], secondary star velocity [vz], 
   orbital inclination [inc], longitude of node [omega] and True anomaly [f]         '''

m_list = get_list(r=0.5, s=0.5, lim=0.5)
x_list = get_list(0.4, 0.1, 5)
vz_list = get_list(-1, s=-1, lim=-1, rel=operator.ge)

inc_list = get_list(r=-((math.pi/2)-math.pi/(inc_count)), s=math.pi/(inc_count), lim=math.pi/2, rel=operator.le)
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
f_list = get_list(r=2*math.pi/f_count, s=2*math.pi/f_count, lim=2*math.pi)

'''Printing the lists lengths to see if there is any problem'''
print(inc_list)
print('length of inclination list= {inclist}'.format(inclist=len(inc_list)))
print('length of Omega list= {omegalist}'.format(omegalist=len(omega_list)))
print('length of True anomaly list= {anolist}'.format(anolist=len(f_list)))

incdifsum = 0
Percentage_array=[]

'''calculating the dimensions of the figure (subplots in height and width)'''
# subnumber = math.ceil(math.sqrt(len(inc_list)))
# print('number of rows and columns of subplots: {sub}'.format(sub=subnumber))

res = {}
with open('/Users/atefeh-behzad/Exoplanet_Simulations/Earth_m_Omega_inc_f/inclination_change_x({min:.2f}-{max:.2f})_{total_count}k.txt'.format(min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as g:
    g.write('inc\tincdifavg\n')
    with open('/Users/atefeh-behzad/Exoplanet_Simulations/Earth_m_Omega_inc_f/Data_x({min:.2f}-{max:.2f})_{total_count}k.txt'.format(min=min(x_list),
                                                                                                                       max=max(x_list),
                                                                                                                       total_count=omega_count*(inc_count)/1000), 'w+') as h:
        h.write('StarM\tStarx\tStarv\tPiinc\tPiOmega\tPif\tPfinc\tPfOmega\tPff\taf\tef\tenergydif%\tangmomdif%\tin?\n')
        for m in m_list:
            print("m= ", m)

            for x in x_list:

                print("x= ", x)

                counter = 0  # this is the counter for calculating percentage

                for v in vz_list:

                    print("vz= ", v)

                    for inc in inc_list:
                        incdifsum = 0
                        inc_name= (inc+math.pi/2)*inc_count/(math.pi)

                        for o in omega_list:
                            print('\nx={x}, Inc={inc:.0f}, Omega: {o:.0f}'.format(x=x, inc=inc_name,o=omega_count*o/(2*math.pi)))
                            print("", end="                         Anomaly # ")
                            for f in f_list:
                                print("{r:.0f}".format(r=f_count*f/(2*math.pi)), end=" ")


                                sim = rebound.Simulation()                              #the simulation starts here

                                sim.add(m=1)                                            #add Sun
                                sim.add(m=3e-6, a=1, Omega=o, inc=inc, f=f)             #add Earth
                                sim.add(m=m, x=x, z=50, vz=v)                           #add the passing star
                                a_c = sim.calculate_orbits()[0]
                                initenergy = sim.calculate_energy()                     #total initial kinetic and potential energy
                                initangmom= sim.calculate_angular_momentum()            #total initial angular momentum in 3 dimensions

                                sim.integrate(round(2 * 50 / (-v), 0))                  #integration

                                tinitangmom= math.sqrt(initangmom[0]**2+initangmom[1]**2+initangmom[2]**2)             #Total size of initial angular momentum

                                try:
                                    angmomdif = math.sqrt(((initangmom[0]-sim.calculate_angular_momentum()[0])**2 + \
                                            (initangmom[1]-sim.calculate_angular_momentum()[1])**2 + \
                                            (initangmom[2]-sim.calculate_angular_momentum()[2])**2)/tinitangmom)        #angular momentum difference devided by initial angular momentum
                                except ZeroDivisionError:
                                    angmomdif = float('inf')


                                a_e = sim.calculate_orbits()[0]                         #calculating orbital elements of the planet after the integration

                                if (abs(a_e.a)*(1+a_e.e))<1.2 and (abs(a_e.a)*(1-a_e.e))>0.8:
                                    Color = '+'
                                else:
                                    Color = ''

                                incdifsum += abs(abs(inc)-abs(a_e.inc))

                                h.write(
                                    '{m:.1f}\t{x:.2f}\t{v}\t{inc:.3f}\t{o:.3f}\t{f:.3f}\t{finc:.3f}\t{fo:.3f}\t{ff:.3f}\t{a:.3f}\t{e:.2g}\t{energychange:.2g}\t\t{angmomdif:.2g}\t\t{Color}\n'.format(m=m, x=x, v=v, inc=inc, o=o, f=f, finc=a_e.inc, fo=a_e.Omega, ff=a_e.f,
                                                                                            a=a_e.a,
                                                                                            e=a_e.e, energychange=(100*abs((initenergy-sim.calculate_energy())/initenergy)), Color=Color, angmomdif=angmomdif))



                        incdifavg= incdifsum/(len(omega_list)*len(f_list))
                        g.write('{inc:.3f}\t{incdifavg:.3f}\n'.format(inc=inc,incdifavg=incdifavg))