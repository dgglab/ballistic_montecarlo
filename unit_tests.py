# you can run these tests with the following command line
# python -m unittest unit_tests.py 
import unittest
import monte_carlo_simulation as mcs
from ballistic_montecarlo.bandstructure.caustic_bandstructure import Bandstructure
from bandstructure.delafossite_bandstructure import delafossite


class TestAll(unittest.TestCase):

  def test_simulation__update_position(self):
    k = delafossite()
    simulation = mcs.Simulation 
    simulation._bandstructure = Bandstructure(k, 0, 10)
    n_f_out, x_out, y_out = simulation._update_position(simulation, [3, 0.6], 8, 8)
    self.assertEqual(n_f_out, (4, 1))
    self.assertEqual(x_out, 8.002412768552901)
    self.assertEqual(y_out, 8.000187869182515)

if __name__ == '__main__':
  unittest.main()
