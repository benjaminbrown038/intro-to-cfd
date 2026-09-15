"""2D MAC projection solver for closed, unit-square cavities.

Nondimensional variables; constant viscosity and thermal diffusivity.
First-order upwind convection and forward Euler time integration.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import factorized


@dataclass
class Config:
    name: str = 'lid_cavity'
    mode: str = 'lid'
    n: int = 32
    reynolds: float = 100.
    prandtl: float = .71
    rayleigh: float = 10000.
    richardson: float = 0.
    end_time: float = 30.
    cfl: float = .4
    dt_max: float = .01
    steady_tolerance: float = 1e-6
    max_steps: int = 200000
    samples: int = 50

    def validate(self):
        if self.mode not in ('lid', 'thermal_lid', 'natural'):
            raise ValueError('mode must be lid, thermal_lid, or natural')
        for key in ('n', 'max_steps', 'samples'):
            v = getattr(self, key)
            if not isinstance(v, int) or isinstance(v, bool) or v < (8 if key == 'n' else 1):
                raise ValueError(f'{key} must be a positive integer (n >= 8)')
        for key in ('reynolds', 'prandtl', 'rayleigh', 'end_time', 'cfl', 'dt_max', 'steady_tolerance'):
            v = getattr(self, key)
            if not np.isfinite(v) or v <= 0:
                raise ValueError(f'{key} must be finite and positive')
        if self.cfl > .5:
            raise ValueError('cfl must be <= 0.5 for this explicit method')
        if not np.isfinite(self.richardson) or self.richardson < 0:
            raise ValueError('richardson must be finite and nonnegative')


class Cavity:
    def __init__(self, config):
        config.validate()
        self.c = config
        self.n, self.h = config.n, 1/config.n
        n = self.n
        self.u = np.zeros((n, n+1))
        self.v = np.zeros((n+1, n))
        self.p = np.zeros((n, n))
        x = (np.arange(n)+.5)/n
        self.theta = np.tile(1-x, (n, 1)) if config.mode != 'lid' else np.full((n,n), .5)
        self.thermal = config.mode != 'lid'
        self.lid = 0. if config.mode == 'natural' else 1.
        if config.mode == 'natural':
            # Free-fall velocity scale sqrt(g beta deltaT L).
            self.nu = np.sqrt(config.prandtl/config.rayleigh)
            self.alpha = 1/np.sqrt(config.rayleigh*config.prandtl)
            self.buoyancy = 1.
        else:
            self.nu = 1/config.reynolds
            self.alpha = self.nu/config.prandtl
            self.buoyancy = config.richardson if self.thermal else 0.
        self.time = 0.
        self.steps = 0
        self.max_energy_defect = 0.
        self.pressure_solve = self._poisson()

    def _poisson(self):
        n, h = self.n, self.h
        a = lil_matrix((n*n, n*n))
        for j in range(n):
            for i in range(n):
                r = j*n+i
                for jj,ii in ((j-1,i),(j+1,i),(j,i-1),(j,i+1)):
                    if 0 <= jj < n and 0 <= ii < n:
                        a[r,jj*n+ii] += 1/h**2
                        a[r,r] -= 1/h**2
        # Pressure reference. Closed-wall flux makes the original RHS compatible.
        a[0,:] = 0
        a[0,0] = 1.
        return factorized(a.tocsc())

    def divergence(self):
        return (np.diff(self.u, axis=1)+np.diff(self.v, axis=0))/self.h

    def project(self, u, v, dt):
        u[:,0] = u[:,-1] = 0
        v[0,:] = v[-1,:] = 0
        rhs = ((np.diff(u,axis=1)+np.diff(v,axis=0))/self.h/dt).ravel()
        rhs[0] = 0.
        p = self.pressure_solve(rhs).reshape((self.n,self.n))
        u[:,1:-1] -= dt*np.diff(p,axis=1)/self.h
        v[1:-1,:] -= dt*np.diff(p,axis=0)/self.h
        return u,v,p-p.mean()

    @staticmethod
    def upwind(speed, backward, forward):
        return np.where(speed >= 0, backward, forward)

    def momentum_rhs(self):
        u,v,h = self.u,self.v,self.h
        # Reflect ghost tangential velocities about the imposed wall velocity.
        ug = np.vstack((-u[0:1], u, 2*self.lid-u[-1:]))
        vg = np.hstack((-v[:,0:1], v, -v[:,-1:]))
        uc,vc = u[:,1:-1],v[1:-1,:]
        vi = .25*(v[:-1,:-1]+v[1:,:-1]+v[:-1,1:]+v[1:,1:])
        ui = .25*(u[:-1,:-1]+u[:-1,1:]+u[1:,:-1]+u[1:,1:])
        ux = self.upwind(uc,(uc-u[:,:-2])/h,(u[:,2:]-uc)/h)
        uy = self.upwind(vi,(uc-ug[:-2,1:-1])/h,(ug[2:,1:-1]-uc)/h)
        vx = self.upwind(ui,(vc-vg[1:-1,:-2])/h,(vg[1:-1,2:]-vc)/h)
        vy = self.upwind(vc,(vc-v[:-2])/h,(v[2:]-vc)/h)
        lapu = (u[:,:-2]+u[:,2:]+ug[:-2,1:-1]+ug[2:,1:-1]-4*uc)/h**2
        lapv = (vg[1:-1,:-2]+vg[1:-1,2:]+v[:-2]+v[2:]-4*vc)/h**2
        ru = -uc*ux-vi*uy+self.nu*lapu
        rv = -ui*vx-vc*vy+self.nu*lapv
        if self.thermal:
            rv += self.buoyancy*(.5*(self.theta[:-1]+self.theta[1:])-.5)
        return ru,rv

    def thermal_rhs(self):
        t,h = self.theta,self.h
        # Total flux = advective - alpha*gradient. All walls impermeable.
        fx = np.zeros_like(self.u)
        fy = np.zeros_like(self.v)
        fx[:,1:-1] = self.u[:,1:-1]*np.where(self.u[:,1:-1]>=0,t[:,:-1],t[:,1:]) - self.alpha*np.diff(t,axis=1)/h
        fy[1:-1] = self.v[1:-1]*np.where(self.v[1:-1]>=0,t[:-1],t[1:]) - self.alpha*np.diff(t,axis=0)/h
        fx[:,0] = 2*self.alpha*(1-t[:,0])/h
        fx[:,-1] = 2*self.alpha*t[:,-1]/h
        rhs = -(np.diff(fx,axis=1)+np.diff(fy,axis=0))/h
        return rhs, float((fx[:,0].sum()-fx[:,-1].sum())*h)

    def timestep(self):
        speed = max(abs(self.u).max(),abs(self.lid))/self.h+abs(self.v).max()/self.h
        diffusivity = max(self.nu,self.alpha if self.thermal else 0.)
        # Boundary ghost diffusion has a larger diagonal; use conservative 8/h².
        rate = speed+8*diffusivity/self.h**2
        return min(self.c.dt_max,self.c.cfl/max(rate,1e-30))

    def step(self, dt):
        old_u,old_v,old_t = self.u.copy(),self.v.copy(),self.theta.copy()
        ru,rv = self.momentum_rhs()
        us,vs = self.u.copy(),self.v.copy()
        us[:,1:-1] += dt*ru
        vs[1:-1] += dt*rv
        self.u,self.v,self.p = self.project(us,vs,dt)
        if self.thermal:
            rhs,net_heat = self.thermal_rhs()
            self.theta += dt*rhs
            defect = abs(float(np.sum(self.theta-old_t)*self.h**2)-dt*net_heat)
            self.max_energy_defect = max(self.max_energy_defect,defect)
        rate = max(abs(self.u-old_u).max(),abs(self.v-old_v).max(),abs(self.theta-old_t).max())/dt
        if not all(np.isfinite(a).all() for a in (self.u,self.v,self.theta,self.p)):
            raise RuntimeError('Nonfinite solution; lower time step or Reynolds/Rayleigh number')
        if self.thermal and (self.theta.min() < -1e-7 or self.theta.max() > 1+1e-7):
            raise RuntimeError('Temperature bounds violated; reduce time step')
        self.time += dt
        self.steps += 1
        return float(rate)

    def fields(self):
        x=(np.arange(self.n)+.5)*self.h
        X,Y=np.meshgrid(x,x)
        u=.5*(self.u[:,:-1]+self.u[:,1:])
        v=.5*(self.v[:-1]+self.v[1:])
        vort=np.gradient(v,self.h,axis=1)-np.gradient(u,self.h,axis=0)
        return dict(x=X,y=Y,u=u,v=v,pressure=self.p.copy(),theta=self.theta.copy(),vorticity=vort)

    def metrics(self, rate=0.):
        f=self.fields()
        return dict(time=self.time,steps=self.steps,change_rate=rate,
                    max_divergence=float(abs(self.divergence()).max()),
                    kinetic_energy=float(.5*np.mean(f['u']**2+f['v']**2)),
                    nusselt_hot=float(np.mean(2*(1-self.theta[:,0])/self.h)) if self.thermal else None,
                    nusselt_cold=float(np.mean(2*self.theta[:,-1]/self.h)) if self.thermal else None,
                    max_energy_step_defect=self.max_energy_defect)

    def run(self):
        history=[]
        frames=[]
        next_sample=0.
        steady_hits=0
        rate=float('inf')
        status='end_time_reached'
        while self.time < self.c.end_time-1e-12:
            if self.steps >= self.c.max_steps:
                status='max_steps_reached'
                break
            dt=min(self.timestep(),self.c.end_time-self.time)
            rate=self.step(dt)
            steady_hits=steady_hits+1 if rate < self.c.steady_tolerance else 0
            if self.time >= next_sample or steady_hits >= 30:
                history.append(self.metrics(rate))
                f=self.fields()
                frames.append(dict(time=self.time,u=f['u'],v=f['v'],theta=f['theta']))
                next_sample=self.time+self.c.end_time/self.c.samples
            if steady_hits >= 30:
                status='steady_tolerance_reached'
                break
        if not history or history[-1]['steps'] != self.steps:
            history.append(self.metrics(rate))
            f=self.fields()
            frames.append(dict(time=self.time,u=f['u'],v=f['v'],theta=f['theta']))
        return dict(config=asdict(self.c),status=status,metrics=self.metrics(rate)),history,frames
