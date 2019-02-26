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
def get_log_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(10**r)
        r += s
    return my_list

'''defining the Maxwell-Boltzman distribution'''
def f_maxwell(v, sigma):
    f_v=0.2437*((v**2)/(sigma*3))*math.exp(-(v**2)*0.22676/(sigma**2))
    return f_v

'''initializing some parameters'''
#####################
x_max=60
x_min0=1.2

v_sigma=120/30
#####################

v_mode=2.1*v_sigma
print('Vmode={v_m:.2f}'.format(v_m=v_mode))
v_min=v_mode-v_sigma
print('V_min={v_m:.2f}'.format(v_m=v_min))
v_max=v_mode+v_sigma
print('V_max={v_m:.2f}'.format(v_m=v_max))
v_step=2*v_sigma/10
print('V_step={v_m:.2f}'.format(v_m=v_step))

vz_list = get_list(v_min, v_step, v_max, rel=operator.le)
print("v:", vz_list)

m_list= get_list(3.1, 0.1, 5.01)
print("m:", m_list)


Percentage_array=[]

# '''calculating the dimensions of the figure (subplots in height and width)'''

for m in m_list:
    res = {}
    with open('/Users/atefeh-behzad/Exoplanet_Simulations/Bulge/Vmode({Vmode:.2f}))_m({m:.2f}).txt'.format(Vmode=v_mode, m=m), 'w+') as g:
        g.write('m\tv\tx\n')
        sum_maxwell = 0  # summing up the maxwell distribution
        sum_sectionmaxwell = 0  # summing up the cross sections on the velocity distribution
        for v in vz_list:
            level = 1
            x_min=x_min0
            while level<6:
                if level==1:
                    omega_count=4
                    inc_count=4
                    f_count=4
                    x_step=0.2
                    print('level 1')
                if level==2:
                    omega_count=5
                    inc_count=5
                    f_count=5
                    x_step=0.02
                    print('level 2')
                if level==3:
                    omega_count=7
                    inc_count=7
                    f_count=7
                    x_step=0.01
                    print('level 3')
                if level==4:
                    omega_count=10
                    inc_count=10
                    f_count=10
                    x_step=0.01
                    print('level 4')
                if level==5:
                    print('level 5')
                    g.write('{m:.2f}\t{v:.2f}\t{x:.2f}\n'.format(m=m, v=v, x=x))
                    # if vx[0,0]==0:
                    #     vx=vxtemp
                    # else:
                    #     np.append(vx,vxtemp)
                    f_v=f_maxwell(v, v_sigma)
                    sum_maxwell+=f_v
                    sum_sectionmaxwell+=(f_v*math.pi*x**2)
                    break

                stepcounter=0
                '''here we make lists of omega,x,vz,inc and m'''
                omega_list = get_list(0, 2 * math.pi / omega_count, 2 * math.pi - math.pi / omega_count,
                                      rel=operator.le)
                inc_list = get_list(r=-((math.pi / 2) - math.pi / (inc_count)), s=math.pi / (inc_count),
                                    lim=math.pi / 2, rel=operator.le)
                f_list = get_list(r=2 * math.pi / f_count - 0.0001, s=2 * math.pi / f_count, lim=2 * math.pi)
                x_list = get_list(x_min, x_step, x_max)
    
                for x in x_list:
    
                    counter = 0                                                         # this is the counter for calculating percentage
                    print('\nLevel={level:.0f}\tv={v:.3f}\tm={m:.2f}\tx={x:.2f}'.format(level=level,v=v, x=x, m=m))
    
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
                                if (abs(a_e.a)*(1+a_e.e))<(1.37) and (abs(a_e.a)*(1-a_e.e))>(0.95):
                                    counter += 1
    
                    Percentage = (round((counter / (inc_count * omega_count * (f_count))) * 100, 2))
                    stepcounter+=1
                    if Percentage >= 95:
                        if stepcounter==1:
                            x_min=x_min-0.05
                            level-=1
                            print("THERE WAS A STEP ERROR")
                        # g.write('{v:.2f}\t{m:.2f}\t{hw:.2f}\t{x:.2f}\n'.format(v=-v, m=m, hw=2*hw, x=x))
                        if level==1:
                            x_min=x-0.8
                        if level==2:
                            x_min=x-0.1
                        if level==3:
                            x_min=x-0.05
                        if level==4:
                            vxtemp=np.array([[v],[x]])
                        level+=1
                        break
        g.write('sum of all maxwell points: {mp:.2f}\nsum of cross sections over maxwell: {cp:.3f}\n'.format(mp=sum_maxwell, cp=sum_sectionmaxwell))
        g.write('mean cross section= {cs:.3f}\n'.format(cs=sum_sectionmaxwell/sum_maxwell))
os.system('say "finished"')