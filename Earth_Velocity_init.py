import rebound
import matplotlib.pyplot as plt
import numpy as np

sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=3e-6, x=1, vy=-1)
sim.t= 0
max_time = 1000000
step_number = 1000
times= np.linspace(1000, max_time, step_number)
with open('/Users/Behzadarbab/Exoplanet_Simulations/Velocity_and_Orbital_test/Data_Velocity_detailed_time{time}.txt'.format(time=max_time), 'w+') as f:
    f.write('a\te\tx0\ty0\tz0\tx1\ty1\tz1\n')
    print('max time is: {time} and we have {steps} steps'.format(time=max_time, steps= step_number))
    print(times)
    for t in times:
        sim.integrate(t)
        print(t)
        a_e = sim.calculate_orbits()[0]
        p0 = sim.particles[0]
        p1 = sim.particles[1]
        f.write('{a}\t{e}\t{x0}\t{y0}\t{z0}\t{x1}\t{y1}\t{z1}\n'.format(a=a_e.a, e=a_e.e, x0=p0.x, y0=p0.y, z0=p0.z,
                                                                        x1=p1.x, y1=p1.y, z1=p1.z))
        # fig = rebound.OrbitPlot(sim, trails=True, slices=True, color=True, periastron=True, unitlabel="[AU]")
        # plt.savefig('/Users/Behzadarbab/Exoplanet_Simulations/Velocity_and_Orbital_test/Velocity_t={time}_of_{times}.png'.format(time=t, times=max_time))
        # plt.close()