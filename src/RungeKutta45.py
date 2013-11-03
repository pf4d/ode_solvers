from scipy.integrate._ode import IntegratorBase
from numpy import *

class RungeKutta45(IntegratorBase):
    runner = True

    def __init__(self,dt=.01,atol=1e-12,rtol=1e-6,S=.98,hmax=10.,hmin=.2):
        self.dt = dt
        self.atol = atol
        self.rtol = rtol
        self.S    = S
        self.hmin = .2
        self.hmax = 10.

    def reset(self,n,has_jac):
        pass

    def run(self,f,jac,y0,t0,t1,f_params,jac_params):
        t = t0
        yo = array(y0)

        c = array([0.,1./5, 3./10, 4./5, 8./9, 1.,1.])
        a = zeros((7,7))
        a[1,0] = 1./5
        a[2,0:2] = [3./40,9./40]
        a[3,0:3] = [44./45,-56./15,32./9]
        a[4,0:4] = [19372./6561,-25360./2187,64448./6561,-212./729]
        a[5,0:5] = [9017./3168,-355./33,46732./5247,49./176,-5103./18656]
        a[6,0:6] = [35./384,0.,500./1113,125./192,-2187./6784,11./84]

        b     = array([35./384, 0., 500./1113, 125./192,-2187./6784,11./84,0.])
        bstar = array([5179./57600,
0.,7571./16695,393./640,-92097./339200,187./2100,1./40])

        k = zeros((7,yo.size))
        
        # Prime k for FSAL
        k[0] = f(t0,y0)
        
        scale = self.atol + self.rtol*abs(yo)
        dnf = sum( (k[0] / scale)**2 )
        dny = sum( (y0 / scale)**2 )
        
        if dnf <= 1e-10 or dny <= 1e-10:
            self.dt = 1e-6
        
        # Perform an explicit Euler step:
        yn = y0 + k[0] * self.dt
        k[1] = f(t0 + self.dt, yn) 
        
        # Estimate the second derivative of the solution:
        der2 = sum( sqrt( abs(k[1] - k[0]) / scale ) ) / self.dt
        
        # step size is computed such that h**5 * max(norm(k[0]),norm(der2)) = 0.01
        der12 = max( abs(der2), sqrt(dnf) )
        if der12 <= 1e-15:
            dtn = max(1e-6, abs(self.dt)*1e-3)
        else:
            dtn = (0.01/der12)**(1/5.)
        self.dt = min( 100.0 * self.dt, min(dtn, self.hmax) )
        
        
        # Integration loop
        while t < t1:
            # Compute ki to include the estimate at y_{n+1}
            # This is for FSAL (first-same-as-last)
            for i in range(1,7):
                yt   = yo + dot(a[i],k) * self.dt
                k[i] = f(t + c[i] * self.dt, yt)

            yn   = yt # 5th order estimate was computed

            # Delta: 5th order minus 4th
            Delta  = dot(b - bstar, k) * self.dt
            
            #Errors
            scale = self.atol + max(abs(yn),abs(yo)) * self.rtol
            err   =  sqrt( sum((Delta / scale)**2) / yn.size)

            # Forward or not, depending of the error values
            if err <= 1.:
                t += self.dt
                self.dt =\
                min(t1-t, self.S * self.dt * (1./err) ** .2, self.hmax*self.dt)
                # note the t1-t above assures final time step hits t1
                yo = yn
                k[0] = k[6] # FSAL assignment

            else:
                self.dt = \
                max(self.S * self.dt * (1./err) ** .2, self.hmin*self.dt)

        if isfinite(yn[-1]): self.success = True
        return yn,t

if RungeKutta45.runner:
    IntegratorBase.integrator_classes.append(RungeKutta45)
