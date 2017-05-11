from scipy.integrate._ode import *
 
import ode_solvers.Euler           as Euler
import ode_solvers.EulerRichardson as EulerRichardson
import ode_solvers.EulerCromer     as EulerCromer
import ode_solvers.Predictor       as Predictor
import ode_solvers.RungeKutta      as RungeKutta
 
from numpy import arange,vstack,array, sqrt
from pylab import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset




# Function defining derivates as a function of the positions and velocities.
# This is the INTERACTION, the inputs to this function are OBSERVABLES.
def sho(t,y,k,m):
  """
  INPUTS:
    t time, only used in non-autonomous systems.
    y[0] the x-position
    y[1] the velocity
    k Spring constant
    m Mass
  OUTPUTS:
    dxdt derivative of the position = velocity
    dvdt acceleration, Newton's second law, here for a spring a = -k/m x
  """
  return array([y[1], -k/m*y[0]])






# Initial Conditions
y0, t0 = [1.,0.], 0.
tf = 200.            # Final time
dt = .2              # Output interval
 
# Model parameters these and the observables define the COMPONENT
k = 1.               # Spring constant
m = 9.               # Mass on the Spring


# Position arrays
yfE = [array(y0)]
yfEC = [array(y0)]
yfER = [array(y0)]
yfRK = [array(y0)]
yfPC = [array(y0)]

# Times to evaluate a solution. 
time = arange(t0,tf,dt)





# CREATE ODE OBJECTS
i = ode(sho)

# Main loops for the integration
# Euler Method:
i.set_integrator('Euler',dt=.05)
i.set_initial_value(y0,t0)
i.set_f_params(k, m)
for t in time:
    i.integrate(i.t+dt)
    yfE.append(i.y)
yfE = array(yfE)

# Euler-Richardson Method:
i.set_integrator('EulerRichardson',dt=.05)
i.set_initial_value(y0,t0)
i.set_f_params(k, m)
for t in time:
    i.integrate(i.t+dt)
    yfER.append(i.y)
yfER = array(yfER)

# Euler-Cromer Method:
i.set_integrator('EulerCromer',dt=.05)
i.set_initial_value(y0,t0)
i.set_f_params(k, m)
for t in time:
    i.integrate(i.t+dt)
    yfEC.append(i.y) 
yfEC = array(yfEC)

# RungeKutta Method:
i.set_integrator('RungeKutta',dt=.05)
i.set_initial_value(y0,t0)
i.set_f_params(k, m)
for t in time:
    i.integrate(i.t+dt)
    yfRK.append(i.y)
yfRK = array(yfRK)

# Predictor-Corrector Method:
i.set_integrator('Predictor',dt=.05)
i.set_initial_value(y0,t0)
i.set_f_params(k, m)
for t in time:
    i.integrate(i.t+dt)
    yfPC.append(i.y)
yfPC = array(yfPC)




# analytic solution
yT = y0[0]*cos(sqrt(k/m)*time)




# The discrepancy  between simulated and analytically calculated time at the
# last simulated position, x[-1]
#disc = time[-1] - sqrt((yfE[-1][0] - y0[0])/(-.5*g))



# Plot the results
purp = '#880cbc'
grun = '#77f343'

fig = figure(figsize=(8,5))
ax = fig.add_subplot(111, ylim=(-2.5,2.5))

#text(0.5, 10, "discrepancy: %f" % (disc))
plot(time, yfE[1:,0],  lw=2.0, color='k',  ls='-',  label='Euler')
#plot(time, yfER[1:,0], lw=2.0, color='b',  ls='-',  label='Euler-Richardson')
plot(time, yfEC[1:,0], lw=2.0, color=grun, ls='-',  label='Euler-Cromer')
#plot(time, yfRK[1:,0], lw=2.0, color=purp, ls='-',  label='Runge-Kutta')
plot(time, yfPC[1:,0], lw=2.0, color=purp, ls='-',  label='Predictor-Corrector')
plot(time, yT,         lw=2.0, color='r',  ls='--', label='Analytical')

# Legend formatting:
leg = legend(loc='upper left')
ltext  = leg.get_texts()
frame  = leg.get_frame()
setp(ltext, fontsize='small')
frame.set_alpha(0.0)
frame.set_facecolor('0.80')

# Main formatting:
xlabel(r'$t$ [s]')
ylabel(r'$x$ [m]')
title('Simple Harmonic Oscillator')
grid()

# Inset axis formatting:
axins = zoomed_inset_axes(ax, 50, loc=3)
x1, x2, y1, y2 = 178.3, 179.7, -1.005, -0.98
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
xticks(visible=False)
yticks(visible=False)
mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5")
plot(time, yfE[1:,0],  lw=2.0, color='k',  ls='-',  label='Euler')
#plot(time, yfER[1:,0], lw=2.0, color='b',  ls='-',  label='Euler-Richardson')
plot(time, yfEC[1:,0], lw=2.0, color=grun, ls='-',  label='Euler-Cromer')
#plot(time, yfRK[1:,0], lw=2.0, color=purp, ls='-',  label='Runge-Kutta')
plot(time, yfPC[1:,0], lw=2.0, color=purp, ls='-',  label='Predictor-Corrector')
plot(time, yT,         lw=2.0, color='r',  ls='--', label='Analytical')

tight_layout()
savefig('sho.pdf')
#show()


