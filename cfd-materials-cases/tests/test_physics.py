"""Independent reference, conservation, limit, and mesh checks."""
import unittest
import numpy as np
from cfd_cases.__main__ import load_cases
from cfd_cases.models import channel, duct, porous, thermal, compliant

C = {c['model']: c for c in load_cases()}


class PhysicsTests(unittest.TestCase):
    def test_channel_second_order(self):
        errors = [channel(dict(C['channel'], n=n))[1]['flow_relative_error'] for n in (10, 20, 40)]
        self.assertAlmostEqual(errors[0]/errors[1], 4., places=7)
        self.assertAlmostEqual(errors[1]/errors[2], 4., places=7)

    def test_inverse_viscosity(self):
        q = channel(C['channel'])[1]['flow_per_depth_m2_s']
        q2 = channel(dict(C['channel'], mu=2*C['channel']['mu']))[1]['flow_per_depth_m2_s']
        self.assertAlmostEqual(q2/q, .5, places=12)

    def test_square_duct_series_and_refinement(self):
        # Exact square-duct series Q=G*a^4/(12 mu)*(1-192/pi^5 sum(tanh(n*pi/2)/n^5)).
        n = np.arange(1, 1000, 2)
        c = dict(C['duct'], shape='rectangle')
        exact = c['pressure_gradient']*c['width']**4/(12*c['mu'])*(1-192/np.pi**5*np.sum(np.tanh(n*np.pi/2)/n**5))
        errors = []
        for size in (16, 32, 64):
            q = duct(dict(c, nx=size, ny=size))[1]['flow_m3_s']
            errors.append(abs(q/exact-1))
        self.assertLess(errors[-1], .002)
        self.assertGreater(errors[0]/errors[1], 3.5)
        self.assertGreater(errors[1]/errors[2], 3.5)

    def test_insert_reduces_flow(self):
        a = duct(dict(C['duct'], shape='rectangle'))[1]['flow_m3_s']
        fields, b = duct(C['duct'])
        self.assertLess(b['flow_m3_s'], a)
        self.assertTrue(np.all(fields['velocity_m_s'][fields['fluid_mask'] == 0] == 0))

    def test_uniform_darcy_reference(self):
        c = dict(C['porous'], pattern='uniform', nx=20, ny=10)
        fields, m = porous(c)
        exact = c['inlet_pressure']*(1-fields['x_m']/c['length'])
        self.assertLess(np.max(abs(fields['pressure_Pa']-exact)), 1e-9)
        q = c['permeability']/c['mu']*c['inlet_pressure']/c['length']*c['height']
        self.assertAlmostEqual(m['inlet_flow_per_depth_m2_s']/q, 1., places=10)

    def test_layered_darcy_resistance(self):
        c = dict(C['porous'], pattern='barrier', nx=50, ny=10)
        _, m = porous(c)
        resistance = c['mu']*c['length']/c['permeability']*(.8+.2/c['contrast'])
        expected = c['height']*c['inlet_pressure']/resistance
        self.assertAlmostEqual(m['inlet_flow_per_depth_m2_s']/expected, 1., places=10)

    def test_porous_mass_conservation(self):
        _, m = porous(C['porous'])
        self.assertLess(m['mass_relative_imbalance'], 1e-9)
        self.assertLess(m['max_divergence_s_inv'], 1e-10)

    def test_thermal_energy_balance_and_bounds(self):
        c = C['thermal']
        f, m = thermal(c)
        self.assertLess(m['energy_relative_imbalance'], 1e-9)
        self.assertGreaterEqual(f['temperature_K'].min(), c['inlet_temperature']-1e-8)
        self.assertLessEqual(f['temperature_K'].max(), c['outer_temperature']+1e-8)

    def test_thermal_equilibrium(self):
        c = dict(C['thermal'], outer_temperature=C['thermal']['inlet_temperature'])
        f, _ = thermal(c)
        self.assertLess(np.max(abs(f['temperature_K']-c['inlet_temperature'])), 1e-8)

    def test_more_conductive_wall(self):
        a = thermal(dict(C['thermal'], k_solid=.2))[1]['outlet_bulk_temperature_K']
        b = thermal(dict(C['thermal'], k_solid=400.))[1]['outlet_bulk_temperature_K']
        self.assertGreater(b, a)

    def test_compliance_against_integrated_polynomial(self):
        c = C['compliant']
        _, m = compliant(c)
        beta = c['radius']/(c['young_modulus']*c['thickness'])
        p0 = c['outlet_pressure']-c['external_pressure']
        p1 = c['inlet_pressure']-c['external_pressure']
        # Integrate the expanded radius^4 exactly, independently of trapezoidal solver.
        integral = sum(a*beta**j*(p1**(j+1)-p0**(j+1))/(j+1) for j,a in enumerate([1,4,6,4,1]))
        exact = np.pi*c['radius']**4*integral/(8*c['mu']*c['length'])
        self.assertLess(abs(m['flow_m3_s']/exact-1), 1e-7)

    def test_rigid_limit_and_pressure_monotonicity(self):
        fields, m = compliant(dict(C['compliant'], young_modulus=1e16))
        self.assertLess(abs(m['flow_m3_s']/m['rigid_flow_m3_s']-1), 1e-10)
        self.assertTrue(np.all(np.diff(fields['pressure_Pa']) < 0))
        self.assertTrue(np.all(np.diff(fields['radius_m']) <= 0))

    def test_invalid_physics_rejected(self):
        with self.assertRaises(ValueError):
            channel(dict(C['channel'], mu=-1))
        with self.assertRaises(ValueError):
            compliant(dict(C['compliant'], young_modulus=1000.))
        with self.assertRaises(ValueError):
            thermal(dict(C['thermal'], n_wall=3))


if __name__ == '__main__':
    unittest.main()
