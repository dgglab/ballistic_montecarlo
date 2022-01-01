# you can run these tests with the following command line
# python -m unittest unit_tests.py 
import unittest
import monte_carlo_simulation as mcs
from ballistic_montecarlo.bandstructure.caustic_bandstructure import Bandstructure
from bandstructure.delafossite_bandstructure import delafossite
from ballistic_montecarlo.geo.caustic_frame import Edge
import geo.caustic_frame as caustic_frame
from copy import deepcopy
import time
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from bandstructure.circular_bandstructure import circular
from matplotlib.animation import FuncAnimation


class TestAll(unittest.TestCase):

  def test_simulation__update_position(self):
    tef_frame = caustic_frame.Frame('geo/tef.dxf')
    k = delafossite()
    simulation = mcs.Simulation(tef_frame, k, 0, 10)
    n_f_out, x_out, y_out = simulation._update_position([3, 0.6], 8, 8)
    self.assertEqual(n_f_out, (4, 1))
    self.assertEqual(x_out, 8.002412768552901)
    self.assertEqual(y_out, 8.000187869182515)

  def test_simulation__specular(self):
    tef_frame = caustic_frame.Frame('geo/tef.dxf')
    k = delafossite()
    simulation = mcs.Simulation(tef_frame, k, 0, 4)
    edge = Edge(2, 4, 7, 3, 0)
    n_f_out = simulation._specular([4, 0.3], edge)
    self.assertEqual(n_f_out, (974, 0.39465996752259236))

  def test_simulation_corner_specular(self):
    tef_frame = caustic_frame.Frame('geo/tef.dxf')
    k = delafossite()
    simulation = mcs.Simulation(tef_frame, k, 0, 4)
    edge0 = Edge(2, 4, 7, 3, 0)
    edge1 = Edge(4, 2, 4, 8, 0)
    n_f_out = simulation._corner_specular([4, 0.3], edge0, edge1, 0)
    self.assertEqual(n_f_out, (985, 0.3081393959289981))

  def test_simulation_get_ts_us(self):
    tef_frame = caustic_frame.Frame('geo/tef.dxf')
    k = delafossite()
    simulation = mcs.Simulation(tef_frame, k, 0, 4)
    ts, us = simulation.get_ts_us(2, 10, 7, 4, 9, 5)
    self.assertEqual(ts, 0.7222222222222222)
    self.assertEqual(us, -1.2222222222222223)

  def test_simulation__get_fermi_intersections(self):
    tef_frame = caustic_frame.Frame('geo/tef.dxf')
    k = delafossite()
    simulation = mcs.Simulation(tef_frame, k, 0, 4)
    edge = Edge(2, 9, 3, 11, 0)
    intersections = simulation._get_fermi_intersections([12, 0.5], edge)
    self.assertEqual(len(intersections), 2)
    self.assertEqual(intersections[0], (12, 0.5000000000003393))
    self.assertEqual(intersections[1], (340, 0.31319998953334305))

if __name__ == '__main__':
  unittest.main()
