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
import os
import numpy


'''we define listing'''
def get_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(r)
        r += s
    return my_list

'''we define listing for vz in 10^x form'''
def get_v_log_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(-10**r)
        r += s
    return my_list

'''here goes the Omega, inclination and true anomaly counts we need'''
omega_count= 20            #number of omegas to try from pi/(omega_count/2) to pi (we dont use pi to 2pi because it is symmetric
inc_count= 20              #number of inclinations to try from -pi/2 to pi/2 (we have 30 steps)
f_count= 20                #number of true anomalies from the step size to 2*pi

'''here we make lists of omega,x,vz,inc and m'''
omega_list = get_list(0, 2* math.pi/omega_count, 2 * math.pi-math.pi/omega_count, rel=operator.le)
x_list = get_list(3.8, 0.05, 100.0)
vz_list = get_v_log_list(0, s=0.04, lim=0, rel=operator.le)
print("v:",vz_list)
inc_list = get_list(r=-((math.pi/2)-math.pi/(inc_count)), s=math.pi/(inc_count), lim=math.pi/2, rel=operator.le)
m_list = get_list(r=1, s=0.2, lim=4)
hwlist = get_list(0.2, s=0.02, lim= 0.2)
f_list = get_list(r=2*math.pi/f_count-0.0001 , s=2*math.pi/f_count, lim=2*math.pi)

'''Printing the lists lengths to see if there is any problem'''
print('length of inclination list= {inclist}'.format(inclist=len(inc_list)))
print('length of Omega list= {omegalist}'.format(omegalist=len(omega_list)))
print('length of True anomaly list= {anolist}'.format(anolist=len(f_list)))

Percentage_array=[]

'''calculating the dimensions of the figure (subplots in height and width)'''
subnumber = math.ceil(math.sqrt(len(inc_list)))
print('number of rows and columns of subplots: {sub}'.format(sub=subnumber))

res = {}
with open('/Users/atefeh-behzad/Exoplanet_Simulations/Half_Chance_Width-Velocity/HZhw({minhw:.2f}-{maxhw:.2f})_x({min:.2f}-{max:.2f})_v_({minv:.2f}-{maxv:.2f})_m{m:.1f}_{total_count}k.txt'.format(minhw=min(hwlist), maxhw=max(hwlist), min=min(x_list),
                                                                                                                       max=max(x_list),minv=min(vz_list), maxv=max(vz_list), m=m,
                                                                                                                       total_count=omega_count*f_count*inc_count/1000), 'w+') as g:
    g.write('HZhw\tv\tx\n')
    for hw in hwlist:

        for m in m_list:
            # print("m= ", m)

            for v in vz_list:

                # flag = 0

                # print("vz= ", v)
                # g.write('v= {v}\n'.format(v=v))

                for x in x_list:
                    # if flag ==0:
                    #     continue
                    # print("x= {x:.2f}".format(x=x))

                    counter = 0  # this is the counter for calculating percentage
                    print('\nhw={hw:.2f}\tv={v:.3f}\tx={x:.2f}'.format(hw=hw, v=v, x=x))

                    for inc in inc_list:

                        inc_name= (inc+math.pi/2)*inc_count/(math.pi)


                        for o in omega_list:
                            # print('\nx={x:.2f}, Inc={inc:.0f}, Omega: {o:.0f}'.format(x=x, inc=inc_name,o=omega_count*o/(2*math.pi)))
                            for f in f_list:

                                sim = rebound.Simulation()                              #the simulation starts here
                                sim.add(m=1)                                            #add Sun
                                sim.add(m=3e-6, a=1, Omega=o, inc=inc, f=f)             #add Earth
                                sim.add(m=m, x=x, z=50, vz=v)                           #add the passing star
                                a_c = sim.calculate_orbits()[0]
                                initenergy = sim.calculate_energy()                     #total initial kinetic and potential energy
                                initangmom= sim.calculate_angular_momentum()            #total initial angular momentum in 3 dimensions
                                sim.integrate(round(2 * 50 / (-v), 0))                  #integration

                                a_e = sim.calculate_orbits()[0]                         #calculating orbital elements of the planet after the integration
                                # star_e = sim.calculate_orbits()[1]


                                '''determine the color of HZ orbits in omega-inc plot'''
                                if (abs(a_e.a)*(1+a_e.e))<(1+hw) and (abs(a_e.a)*(1-a_e.e))>(1-hw):
                                    counter += 1

                    Percentage = (round((counter / (inc_count * omega_count * (f_count))) * 100, 2))
                    if Percentage >= 50:
                        g.write('{HZhwidth:.2f}\t{v:.2f}\t{x:.2f}\n'.format(HZhwidth=hw, v=-v, x=x))
                        # g.write('star perigee = {star_p}'.format(star_p=(abs(a_e.a)*(1-a_e.e))))
                        break

os.system('say "A simulation has finished"')