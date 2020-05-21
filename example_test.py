import unittest

import numpy as np

import geo.caustic_frame as caustic_frame
from bandstructure.delafossite_bandstructure import delafossite
import monte_carlo_simulation as mcs


class TestMCS(unittest.TestCase):
    def test_mcs(self):
        k = delafossite()
        bar_frame = caustic_frame.Frame('geo/bar.dxf')
        fields = np.linspace(10, 11, 2)

        for field in fields:
            np.random.seed(42)
            bar_sim = mcs.Simulation(bar_frame, k, 0, field)
            edge_to_collisions, trajectories = bar_sim.run_simulation(2)
            print(f'{field:.1f} T: {list(map(len, trajectories))}')
            self.assertEqual(len(trajectories), 2)


if __name__ == '__main__':
    unittest.main()
