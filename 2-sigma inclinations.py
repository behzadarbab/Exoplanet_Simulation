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
x_step = 0.01
omega_count = 10
inc_count = 1
f_count = 10

#####################
x_max = 7.001
x_min = 3
v = 1
m=1
#####################

res = {}
with open('/Users/atefeh-behzad/Exoplanet_Simulations/2-Sigma/inclination_v({v:.2f})_m({m:.2f})_x({xmin:.2f}_{xmax:.2f}).txt'.format(v=v, m=m, xmin=x_min, xmax=x_max) , 'w+') as g:
    g.write('inc\tx\tpercent\n')

    '''here we make lists of omega,x,vz,inc and m'''
    omega_list = get_list(0, 2 * math.pi / omega_count, 2 * math.pi - math.pi / omega_count,
                          rel=operator.le)
    inc_list = get_list(r=-((math.pi / 2) - math.pi / (inc_count)), s=math.pi / (inc_count),
                        lim=math.pi / 2, rel=operator.le)
    f_list = get_list(r=2 * math.pi / f_count - 0.0001, s=2 * math.pi / f_count, lim=2 * math.pi)
    x_list = get_list(x_min, x_step, x_max)

    for x in x_list:

        print('\nv={v:.2f}\tx={x:.3f}'.format(v=v, x=x))

        for inc in inc_list:
            inc=0.01
            counter = 0  # this is the counter for calculating percentage
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
                    if (abs(a_e.a) * (1 + a_e.e)) < (1.37) and (abs(a_e.a) * (1 - a_e.e)) > (0.95):
                        counter += 1
            Percentage = (round((counter / (omega_count * f_count)) * 100, 4))
            g.write('{inc:.2f}\t{x:.2f}\t{percent:.4f}\n'.format(inc=inc, x=x, percent=Percentage))
os.system('say "finished"')