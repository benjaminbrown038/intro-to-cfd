import unittest
import numpy as np
from cavity_cases.mesh_study import field_changes,relative_change,validate_levels,assess


def manufactured(n):
    x,y=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    return dict(x=x,y=y,u=2*x-y,v=x+3*y,theta=1-x)


class MeshStudyTests(unittest.TestCase):
    def test_linear_fields_interpolate_exactly(self):
        r=field_changes(manufactured(16),manufactured(32))
        self.assertLess(r['velocity_sampled_l2_change_percent'],1e-12)
        self.assertLess(r['temperature_rms_change_percent_of_deltaT'],1e-12)

    def test_known_velocity_change(self):
        coarse=manufactured(16)
        coarse['u']*=1.1; coarse['v']*=1.1
        r=field_changes(coarse,manufactured(48))
        self.assertAlmostEqual(r['velocity_sampled_l2_change_percent'],10.,places=10)

    def test_grids_and_zero_reference(self):
        self.assertEqual(validate_levels([48,16,32]),[16,32,48])
        for levels in ([16,32],[16,16,32],[4,16,32]):
            with self.assertRaises(ValueError): validate_levels(levels)
        self.assertEqual(relative_change(0,0),0.)
        self.assertIsNone(relative_change(1,0))
        self.assertEqual(relative_change(9,10),10.)

    def test_nonsteady_flag(self):
        row=dict(status='end_time_reached',kinetic_energy_change_percent=.1,velocity_sampled_l2_change_percent=.1)
        r=assess([row],1.,False)
        self.assertFalse(r['all_runs_steady'])
        self.assertIn('inconclusive',r['assessment'])

if __name__=='__main__': unittest.main()
