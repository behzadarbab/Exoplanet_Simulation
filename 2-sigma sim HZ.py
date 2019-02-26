"""                                             Briefing about the simulation
This program simulates a three body system; Sun, Earth, a passing star. The passing star starts from 50 AU above the
Solar system plane. We throw the star with a speed (or we can change the speed or give a list of speeds), then after a
time that it reaches around 50 AU below the solar plane, we stop the single simulation and check if the planet orbit
reaches beyond the HZ aroound sun (0.8AU-1.2AU). Then we do this simulation for different inclinations and directions
of the planet orbit and we change the place of the planet on the orbit. After this, we can make a 3D plot (a grand plot
with many subplots in 2D) and if in a simulation, after the flyby, the planet stays in HZ, we plot it in black, otherwise
in gray."""

"""
this is a different approach using a 3-level search
"""


import rebound
import math
import operator
import os
import numpy as np


'''we define listing'''


def get_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(r)
        r += s
    return my_list


'''we define listing for vz in 10**x form'''


def get_ln_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(np.e**r)
        r += s
    return my_list


'''initializing some parameters'''
x = 1
x_step = 0.1

#####################
x_max = 100
x_min0 = 2
v = 3.0
m=1
#####################

hz_list = get_ln_list(-4.0, s=0.2, lim=-0.19, rel=operator.le)
print("hz:", hz_list)

Percentage_array = []

res = {}
with open('/Users/atefeh-behzad/Exoplanet_Simulations/2-Sigma/HZ_hz({minhz:.3f}-{maxhz:.3f})_v({v:.2f}_m({m:.2f}).txt'.format(minhz=min(hz_list), maxhz=max(hz_list), v=v, m=m) , 'w+') as g:
    g.write('hz\tLog(v)\tLog(x)\n')
    for hz in hz_list:
        level = 0
        x_min = x_min0
        while level < 6:
            if level == 0:
                omega_count = 3
                inc_count = 3
                f_count = 3
                x_step = 1
                print('level 0')
            if level == 1:
                omega_count = 4
                inc_count = 4
                f_count = 4
                x_step = 0.2
                print('level 1........................................')
            if level == 2:
                omega_count = 5
                inc_count = 5
                f_count = 5
                x_step = 0.1
                print('level 2........................................')
            if level == 3:
                omega_count = 8
                inc_count = 8
                f_count = 8
                x_step = 0.04
                print('level 3........................................')
            if level == 4:
                omega_count = 10
                inc_count = 10
                f_count = 10
                x_step = 0.01
                print('level 4........................................')
            if level == 5:
                print('level 5........................................')
                g.write('{hz:.4f}\t{v:.2f}\t{Lnx:.4f}\n'.format(hz=hz, v=v, Lnx=np.log(x)))
                break
            '''here we make lists of omega,x,vz,inc and m'''
            omega_list = get_list(0, 2 * math.pi / omega_count, 2 * math.pi - math.pi / omega_count,
                                  rel=operator.le)
            inc_list = get_list(r=-((math.pi / 2) - math.pi / (inc_count)), s=math.pi / (inc_count),
                                lim=math.pi / 2, rel=operator.le)
            f_list = get_list(r=2 * math.pi / f_count - 0.0001, s=2 * math.pi / f_count, lim=2 * math.pi)
            x_list = get_list(x_min, x_step, x_max)

            for x in x_list:

                counter = 0                                                         # this is the counter for calculating percentage
                print('\nhz={hz:.3f}\tv={v:.2f}\tx={x:.3f}'.format(v=v, x=x, hz=hz))

                for inc in inc_list:

                    inc_name= (inc+math.pi/2)*inc_count/(math.pi)

                    for o in omega_list:

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

                            '''determine the color of HZ orbits in omega-inc plot'''
                            if (abs(a_e.a)*(1+a_e.e))<(1+hz) and (abs(a_e.a)*(1-a_e.e))>(1-hz):
                                counter += 1

                Percentage = (round((counter / (inc_count * omega_count * f_count)) * 100, 2))
                if Percentage >= 95:
                    if level == 0:
                        x_min = x-2
                    if level == 1:
                        x_min = x-1.0
                    if level == 2:
                        x_min = x-0.25
                    if level == 3:
                        x_min = x-0.1
                    level += 1
                    break
os.system('say "finished"')