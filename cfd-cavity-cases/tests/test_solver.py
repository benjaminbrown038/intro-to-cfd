import unittest
import numpy as np
from cavity_cases.solver import Config,Cavity


class SolverTests(unittest.TestCase):
    def test_projection_mass_and_boundary(self):
        s=Cavity(Config(n=16))
        rng=np.random.default_rng(12)
        s.u,s.v,s.p=s.project(rng.normal(size=s.u.shape),rng.normal(size=s.v.shape),.001)
        self.assertLess(abs(s.divergence()).max(),1e-9)
        self.assertTrue(np.all(s.u[:,[0,-1]]==0))
        self.assertTrue(np.all(s.v[[0,-1],:]==0))

    def test_projection_removes_gradient(self):
        s=Cavity(Config(n=16))
        rng=np.random.default_rng(8)
        phi=rng.normal(size=s.p.shape)
        u,v=s.u.copy(),s.v.copy()
        u[:,1:-1]=np.diff(phi,axis=1)/s.h
        v[1:-1]=np.diff(phi,axis=0)/s.h
        u,v,p=s.project(u,v,.1)
        self.assertLess(max(abs(u).max(),abs(v).max()),1e-9)

    def test_motionless_equilibrium(self):
        s=Cavity(Config(n=12))
        s.lid=0
        for _ in range(10):
            s.step(s.timestep())
        self.assertEqual(np.max(abs(s.u))+np.max(abs(s.v)),0.)

    def test_linear_conduction_and_unit_nusselt(self):
        s=Cavity(Config(mode='thermal_lid',n=16,richardson=0))
        s.lid=0
        initial=s.theta.copy()
        for _ in range(20):
            s.step(s.timestep())
        self.assertLess(abs(s.theta-initial).max(),1e-12)
        self.assertAlmostEqual(s.metrics()['nusselt_hot'],1.,places=12)

    def test_conduction_second_order_eigenmode(self):
        errors=[]
        for n in (16,32,64):
            s=Cavity(Config(mode='thermal_lid',n=n))
            x,y=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
            mode=.1*np.sin(np.pi*x)*np.cos(np.pi*y)
            s.theta=1-x+mode
            numerical,_=s.thermal_rhs()
            exact=-2*np.pi**2*s.alpha*mode
            errors.append(np.linalg.norm(numerical-exact)/np.linalg.norm(exact))
        self.assertGreater(errors[0]/errors[1],3.9)
        self.assertGreater(errors[1]/errors[2],3.9)

    def test_thermal_flux_balance(self):
        s=Cavity(Config(mode='thermal_lid',n=16))
        rng=np.random.default_rng(5)
        s.theta=rng.uniform(.1,.9,s.theta.shape)
        for _ in range(50):
            s.step(s.timestep())
        self.assertLess(s.max_energy_defect,1e-12)
        self.assertGreaterEqual(s.theta.min(),0)
        self.assertLessEqual(s.theta.max(),1)

    def test_hot_side_rises(self):
        s=Cavity(Config(mode='natural',n=16))
        for _ in range(100):
            s.step(s.timestep())
        self.assertGreater(s.v[8,2],0)
        self.assertLess(s.v[8,-3],0)
        self.assertLess(abs(s.divergence()).max(),1e-9)

    def test_invalid_inputs(self):
        for kwargs in ({'n':0},{'cfl':1},{'prandtl':-1},{'mode':'unknown'}):
            with self.assertRaises(ValueError):
                Cavity(Config(**kwargs))

    def test_incomplete_status(self):
        s=Cavity(Config(n=8,max_steps=1))
        r,history,frames=s.run()
        self.assertEqual(r['status'],'max_steps_reached')
        self.assertEqual(history[-1]['steps'],1)

if __name__=='__main__':
    unittest.main()
