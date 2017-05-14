import rebound
import math
import operator

def get_list(r, s, lim, rel=operator.le):
    my_list = []
    while rel(r, lim):
        my_list.append(r)
        r += s
    return my_list

omega_list = get_list(math.pi/10,math.pi/10, 2*math.pi)
x_list = get_list(2, 1, 30)
vz_list = get_list(-1, s=-0.5, lim=-1, rel=operator.ge)
inc_list = get_list(r=0, s=math.pi/20, lim=0)
m_list = get_list(r=0.5, s=0.5, lim=10)

res = {}
with open('/Users/Behzadarbab/Data_Earth_Simple.txt', 'w+') as f:
    for m in m_list:
        print ("m= ", m)
        for x in x_list:
            print("x= ", x)
            for v in vz_list:
                for inc in inc_list:
                    for o in omega_list:
                        sim = rebound.Simulation()
                        sim.add(m=1)
                        sim.add(m=3e-6, a=1, Omega=o, inc=inc)
                        sim.add(m=m, x=x, z=150, vz=v)
                        sim.integrate(2 * 150/(-v))
                        a_e = sim.calculate_orbits()[0]

                        res[(m, x, v, inc, o)] = (a_e.a, a_e.e)
                        f.write('{m}, {x}, {v}, {inc:.3f}, {o:.3f},  {a:.4f}, {e:.4f}\n'.format(m=m, x=x, v=v, inc=inc, o=o, a=a_e.a,
                                                                              e=a_e.e))

