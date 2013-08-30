from scipy.integrate.ode import IntegratorBase
from numpy import array,arange,isfinite,ceil
from pylab import linspace

class RungeKutta(IntegratorBase):
    runner = True

    def __init__(self,dt=.01):
        self.dt = dt

    def reset(self,n,has_jac):
        pass

    def run(self,f,jac,y0,t0,t1,f_params,jac_params):
      # this method is called to integrate from t=t0 to t=t1
      # with initial condition y0. f and jac are user-supplied functions
      # that define the problem. f_params,jac_params are additional
      # arguments to these functions.

     yo = array(y0) # Initial condition
     t = t0
     dt = self.dt

     while t < t1:
     
       # For last value of dt:
       if t + dt > t1:
         dt = t1 - t
         
       k1 = f(t, yo, *f_params) * dt
       k2 = f(t + dt/2., yo + k1/2., *f_params) * dt
       k3 = f(t + dt/2., yo + k2/2., *f_params) * dt
       k4 = f(t + dt, yo + k3, *f_params) * dt
       yn = yo.copy() + (1./6.) * (k1 + 2*k2 + 2*k3 + k4)
       yo = yn.copy()
       
       t += dt

     if isfinite(yn[-1]): self.success = 0 # Check for success
     return yn,t

if RungeKutta.runner:
    IntegratorBase.integrator_classes.append(RungeKutta)
