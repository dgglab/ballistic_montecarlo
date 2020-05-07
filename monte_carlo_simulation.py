import numpy as np
from bandstructure.caustic_bandstructure import Bandstructure
from shapely.geometry import LineString


class Simulation:
    def __init__(self, frame, k, phi, B, n_inject):
        self.frame = frame
        self.bandstructure = Bandstructure(k, phi, B)
        self.B = B
        self.phi = phi
        self.n_inject = n_inject

        self.calculate_injection_probs()

    def calculate_injection_probs(self):
        for edge in self.frame.edges:
            in_prob, cum_prob = self.bandstructure.calculate_injection_prob(
                edge.normal)
            edge.set_in_prob(in_prob, cum_prob)

    def run_simulation(self):
        x_store = []
        y_store = []
        n_store = []

        for _ in range(self.n_inject):
            x_cur_store = []
            y_cur_store = []
            n_cur_store = []

            (x, y), injecting_edge = self.frame.get_inject_position(1)
            n = injecting_edge.get_injection_index()

            x_cur_store.append(x)
            y_cur_store.append(y)
            n_cur_store.append(n)

            active_trajectory = True
            while active_trajectory:
                [x_new, y_new], n_new = self.update_position([x, y], n)

                # Can do floating point calc to save time here
                intersection = False
                intersection_edges = []
                for edge in self.frame.edges:
                    step = LineString([(x, y), (x_new, y_new)])
                    edge_int = step.intersects(edge.linestring)
                    if edge_int:
                        intersection = True
                        intersection_edges.append(edge)

                # this code will not work if there are multiple intersections beyond the start point
                if intersection:
                    for edge in intersection_edges:
                        x_int, y_int = step.intersection(
                            edge.linestring).coords.xy
                        if x_int[0] != x or y_int[0] != y:
                            x_new = x_int[0]
                            y_new = y_int[0]

                            if edge.layer == 0:
                                n_new = self.scatter(edge)
                            if edge.layer == 1:
                                n_new = self.scatter(edge)
                            if edge.layer == 2:
                                active_trajectory = False

                x = x_new
                y = y_new
                n = n_new
                x_cur_store.append(x)
                y_cur_store.append(y)
                n_cur_store.append(n)

            x_store.append(x_cur_store)
            y_store.append(y_cur_store)
            n_store.append(n_cur_store)

        return x_store, y_store, n_store

    def update_position(self, r, n):
        r = r + self.bandstructure.dr[:, n]
        n = (n + 1) % (np.shape(self.bandstructure.dr)[1])
        return r, n

    def scatter(self, edge):
        return edge.get_injection_index()
