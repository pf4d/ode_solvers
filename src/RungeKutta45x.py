from scipy.integrate.ode import IntegratorBase
from numpy import *

class RungeKutta45(IntegratorBase):
    runner = True

    def __init__(self,dt=.01):
      self.dt = dt
      self.a = array \
      (
        [
        [0., 0., 0., 0., 0., 0., 0.], 
        [1/5., 0., 0., 0., 0., 0., 0.], 
        [3/40., 9/40., 0., 0., 0., 0., 0.], 
        [44/45., -56/15., 32/9., 0., 0., 0., 0.], 
        [19372/6561., -25360/2187., 64448/6561., -212/729., 0., 0., 0.], 
        [9017/3168., -355/33., 46732/5247., 49/176., -5103/18656., 0., 0.], 
        [35/384., 0., 500/1113., 125/192., -2187/6784., 11/84., 0.],
        [5179/57600., 0., 7571/16695., 393/640.,-92097/339200.,187/2100.,0.]
        ]
      )
      self.b = self.a[6,:]
      self.b_star = self.a[7,:]
      self.c = array \
      (
        [
        [0.], 
        [1/5.], 
        [3/10.], 
        [4/5.], 
        [8/9.], 
        [1.], 
        [1.]
        ]
      )


    def reset(self,n,has_jac):
      pass
    
    
    # this method is called to integrate from t=t0 to t=t1
    # with initial condition y0. f and jac are user-supplied functions
    # that define the problem. f_params,jac_params are additional
    # arguments to these functions.
    def run(self,f,jac,y0,t0,t1,f_params,jac_params):
    
      # Initial conditions
      y0 = array(y0)
      atol = 1e-6
      rtol = 1e-6
      t = t0
      err = 2
      facold = 1e-4
      
      a = self.a
      b = self.b
      c = self.c
      b_star = self.b_star
      dt = self.dt
      k = zeros( (7,len(y0)) )
      

      # Integration loop
      while t < t1:
        
        # For last value of dt, return the fourth-order estimate:
        if t + dt > t1:     
          dt = t1 - t
          for i in range(7):
            k[i] = f(t + c[i]*dt, y0 + dot(a[i,:], k) * dt, *f_params)
          yn = y0 + dot(b, k)
        
        # Otherwise, continue integrating until the error between 
        # the fifth-order and fourth-order formulas for a given dt is 1:
        #   yn      = fourth-order
        #   yn_star = fifth-order
        #   err     = error between formulas
        #   dtn     = new time step   
        else:
        
          while err > 1:
            
            for i in range(7):
              k[i] = f(t + c[i]*dt, y0 + dot(a[i,:], k) * dt, *f_params)
            
            yn = dot(b, k)
            yn_star = dot(b_star, k)
            
            delta = dt * (yn - yn_star)
            
            y = maximum( abs(y0), abs(yn) )
            scale = atol + y*rtol
            
            err = sqrt( sum(delta/scale)**2 / len(yn) )
            
            fac11 = err**(0.2-0.04*0.75)
            fac = fac11 / facold**0.04
            fac = max(1/10., min(1/0.2, fac/0.9))
            dtn = dt / fac
            facold = max(err, 1e-4)
            
            if dtn >= t1-t :
              dtn = t1-t
            elif dtn <= dt * 1e-10:
              dtn = dt * 1e-10
            self.dt = dtn
          
          err = 2
          t += dt
          y0 = y0 + yn
      
      if isfinite(yn[-1]): self.success = 0 # Check for success
      return yn,t

if RungeKutta45.runner:
    IntegratorBase.integrator_classes.append(RungeKutta45)
