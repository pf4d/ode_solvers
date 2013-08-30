#!/usr/bin/env python
from scipy.integrate.ode import *
import time

import Euler as Euler
import EulerRichardson as EulerRichardson
import RungeKutta as RungeKutta
import RungeKutta45 as RungeKutta45

from pylab import *
from scipy.special import erf

def rhs(t, x):
  return 2/sqrt(math.pi) * math.exp(-t**2)

integrators = ['Euler', 'EulerRichardson', 
               'RungeKutta', 'RungeKutta45', 'dopri5']
lineStyle   = ['go-', 'ro-', 'bo-', 'yo-', 'k-']

dt  = 0.01
dts = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
c   = 0

intTime = []  # computation time array
for i in range(len(integrators)): intTime.append([])
for i in integrators:
  errs = []
  for t in dts:
    o = ode(rhs)
    if i == 'dopri5':
      o.set_integrator(i)
    else:
      o.set_integrator(i, dt=t)
    o.set_initial_value(0, 0)
    t1 = time.time()
    o.integrate(1.0)
    t2 = time.time()
    intTime[c].append(t2 - t1)
    print i + ' time: %f\n' % (t2 - t1)
    diff = abs(erf(1) - o.y[0])
    errs.append(diff)
  plot(dts, errs, lineStyle[c], label=i)
  c += 1

# plot the errors :
grid()
legend(loc='lower right')
xlabel('Times (s)')
ylabel('Log Errors')
semilogy()
show()

# plot the computation time :
figure()
for i, int in zip(range(len(integrators)), integrators):
  plot(dts, intTime[i], 'x-', label=int)
grid()
title('Computation Time')
xlabel('Time-Step')
ylabel('Time to Compute (s)')
legend()
show()


