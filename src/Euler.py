from scipy.integrate._ode import IntegratorBase
from numpy import array, isfinite

class Euler(IntegratorBase):
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
        
        yn = yo + f(t,yo,*f_params) * dt
        yo = yn.copy()
        t += dt

      if isfinite(yn[-1]): self.success = True # Check for success
      return yn,t

if Euler.runner:
    IntegratorBase.integrator_classes.append(Euler)
