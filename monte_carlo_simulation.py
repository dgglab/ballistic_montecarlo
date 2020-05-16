import time
import numpy as np
from bandstructure.caustic_bandstructure import Bandstructure
from geo.caustic_frame import Edge
from shapely.geometry import LineString
from shapely.geometry import Point


class Simulation:
    def __init__(self, frame, k, phi, B, n_inject, p_scatter=1, p_ohmic_absorb=1):
        self.frame = frame
        self.bandstructure = Bandstructure(k, phi, B)
        self.B = B
        self.phi = phi
        self.n_inject = n_inject
        self.p_scatter = p_scatter
        self.p_ohmic_absorb = p_ohmic_absorb

        self.set_seed(time.time())

        self.calculate_injection_probs()

    def set_seed(self, seed):
        # Will want to use generators when parallelizing
        np.random.seed(int(seed))
        self.rng_seed = seed

    def calculate_injection_probs(self):
        for edge in self.frame.edges:
            in_prob, cum_prob = self.bandstructure.calculate_injection_prob(
                edge.normal_angle)
            edge.set_in_prob(in_prob, cum_prob)

    def run_simulation(self, debug=False):
        trajectories = []
        self.debug = debug

        for _ in range(self.n_inject):
            trajectory = []

            (x, y), injecting_edge = self.frame.get_inject_position(1)
            n_f = injecting_edge.get_injection_index()

            step_params = [(n_f, x, y)]
            trajectory.append(step_params[0])

            active_trajectory = True
            while active_trajectory:
                n_f = step_params[-1][0]
                x = step_params[-1][1]
                y = step_params[-1][2]
                step_params, active_trajectory = self.step_position(n_f, x, y)
                trajectory.extend(step_params)

            trajectories.append(trajectory)

        return trajectories

    def step_position(self, n_f, x, y):
        step_params = []
        active_trajectory = True
        if self.debug and not self.frame.body.intersects(Point(x, y)):
            print('Previous step stepped out of bounds')
            print(n_f, x, y)
            return [(n_f, x, y)], False

        n_f_new, x_new, y_new = self.update_position(n_f, x, y)

        step_coords = ([(x, y), (x_new, y_new)])
        intersections = self.get_sorted_intersections(step_coords)

        if len(intersections) == 0:
            pass
            step_params.append((n_f_new, x_new, y_new))
            # Single edge intersection
        elif len(intersections) == 1 or intersections[0][3] != intersections[1][3]:
            edge, x_int, y_int, _ = intersections[0]

            if edge.layer == 0:  # Device edge
                edge.num_collisions += 1
                r = np.random.rand()
                if r < self.p_scatter:
                    n_f_new = self.scatter(edge)
                    step_params.append((n_f_new, x_int, y_int))
                else:
                    # replace with specular reflection
                    n_f_new = self.scatter(edge)
                    step_params.append((n_f_new, x_int, y_int))

            elif edge.layer == 2:  # Grounded ohmic
                r = np.random.rand()
                if r < self.p_ohmic_absorb:  # absorbed
                    edge.num_collisions += 1
                    step_params.append((n_f_new, x_int, y_int))
                    active_trajectory = False
                else:
                    r = np.random.rand()
                    if r < self.p_scatter:
                        n_f_new = self.scatter(edge)
                        step_params.append((n_f_new, x_int, y_int))
                    else:
                        # replace with specular reflection
                        n_f_new = self.scatter(edge)
                        step_params.append((n_f_new, x_int, y_int))

            else:  # Generic ohmic
                r = np.random.rand()
                if r < self.p_ohmic_absorb:  # absorb and reemit
                    edge.num_collisions += 1
                    step_params.append((n_f_new, x_int, y_int))
                    (x_new, y_new), reinjecting_edge = self.frame.get_inject_position(
                        edge.layer)
                    n_f_new = reinjecting_edge.get_injection_index()
                    step_params.append((n_f_new, x_new, y_new))
                else:
                    r = np.random.rand()
                    if r < self.p_scatter:
                        n_f_new = self.scatter(edge)
                        step_params.append((n_f_new, x_int, y_int))
                    else:
                        # replace with specular reflection
                        n_f_new = self.scatter(edge)
                        step_params.append((n_f_new, x_int, y_int))

        # Corner intersection
        else:
            edge_0, x_new, y_new, _ = intersections[0]
            edge_1, _, _, _ = intersections[1]

            edges = [edge_0, edge_1]
            index = np.argmax([edge_0.layer, edge_1.layer])
            layer = edges[index].layer
            edge_for_count = edges[index]  # to avoid double counting

            if layer == 0:  # Device edge
                edge_for_count.num_collisions += 1
                r = np.random.rand()
                if r < self.p_scatter:
                    n_f_new = self.corner_scatter(edge_0, edge_1)
                    step_params.append((n_f_new, x_int, y_int))
                else:
                    # replace with specular reflection
                    n_f_new = self.corner_scatter(edge_0, edge_1)
                    step_params.append((n_f_new, x_int, y_int))

            elif layer == 2:  # Grounded ohmic
                r = np.random.rand()
                if r < self.p_ohmic_absorb:
                    edge_for_count.num_collisions += 1
                    step_params.append((n_f_new, x_int, y_int))
                    active_trajectory = False
                else:
                    r = np.random.rand()
                    if r < self.p_scatter:
                        n_f_new = self.corner_scatter(edge_0, edge_1)
                        step_params.append((n_f_new, x_int, y_int))
                    else:
                        # replace with specular reflection
                        n_f_new = self.corner_scatter(edge_0, edge_1)
                        step_params.append((n_f_new, x_int, y_int))

            else:  # Generic ohmic
                r = np.random.rand()
                if r < self.p_ohmic_absorb:  # absorb and reemit
                    edge_for_count.num_collisions += 1
                    step_params.append((n_f_new, x_int, y_int))
                    (x_new, y_new), reinjecting_edge = self.frame.get_inject_position(layer)
                    n_f_new = reinjecting_edge.get_injection_index()
                    step_params.append((n_f_new, x_new, y_new))
                else:
                    r = np.random.rand()
                    if r < self.p_scatter:
                        n_f_new = self.corner_scatter(edge_0, edge_1)
                        step_params.append((n_f_new, x_int, y_int))
                    else:
                        # replace with specular reflection
                        n_f_new = self.corner_scatter(edge_0, edge_1)
                        step_params.append((n_f_new, x_int, y_int))

        return step_params, active_trajectory

    def update_position(self, n_f_in, x_in, y_in):
        [x_out, y_out] = [x_in, y_in] + n_f_in[1] * \
            self.bandstructure.dr[:, n_f_in[0]]
        n_f_out = ((n_f_in[0] + 1) % (np.shape(self.bandstructure.dr)[1]), 1)
        return n_f_out, x_out, y_out

    def scatter(self, edge):
        return edge.get_injection_index()

    def corner_scatter(self, edge_0, edge_1):
        # Convolve the probility distributions
        in_prob = (edge_0.in_prob*edge_1.in_prob) / \
            np.sum(edge_0.in_prob*edge_1.in_prob)
        cum_prob = np.cumsum(in_prob)
        return Edge.compute_injection_index(cum_prob)

    def get_sorted_intersections(self, step_coords):
        '''
        Calls get intersctions and sorts the returned intersctions by their distance
        '''
        # TODO coords is slow!!!
        x = step_coords[0][0]
        y = step_coords[0][1]
        intersections = self.get_intersections(step_coords)

        def compute_s(intersection):
            edge, x_int, y_int = intersection
            S = np.sqrt((x_int - x)**2 + (y_int - y)**2)
            return (edge, x_int, y_int, S)

        intersections_with_S = map(compute_s, intersections)
        if len(intersections) < 2:
            return list(intersections_with_S)
        else:
            return sorted(intersections_with_S, key=lambda intersections: intersections[3])

    def get_intersections(self, step_coords, bias=True):
        '''
        Given a step, check all the edges in the frame to see if the step crosses
        Return True if it does and a list of crossed edges
        '''
        x = step_coords[0][0]
        y = step_coords[0][1]
        x_new = step_coords[1][0]
        y_new = step_coords[1][1]

        x_del = x_new - x
        y_del = y_new - y
        x01 = -x_del
        y01 = -y_del
        x02 = x - self.frame.px0
        y02 = y - self.frame.py0

        intersections = []

        ts = (x02*self.frame.y23 - y02*self.frame.x23) / \
             (x01*self.frame.y23 - y01*self.frame.x23)
        us = -(x01*y02 - y01*x02) / (x01*self.frame.y23 - y01*self.frame.x23)

        for i, (t, u) in enumerate(zip(ts, us)):
            if 0 <= t and t <= 1 and 0 <= u and u <= 1:
                edge = self.frame.edges[i]
                if x_del * edge.normal[0] + y_del * edge.normal[1] < 0:
                    x_int = edge.xs[0] + u*(edge.xs[1] - edge.xs[0])
                    y_int = edge.ys[0] + u*(edge.ys[1] - edge.ys[0])
                    if x_int != x or y_int != y:
                        if bias:
                            bias_vector = 1E-10 * \
                                np.array([(x_new-x), (y_new-y)]) / \
                                np.sqrt((x_new-x)**2 + (y_new-y)**2)
                            intersections.append(
                                (edge, x_int - bias_vector[0], y_int - bias_vector[1]))

                        else:
                            intersections.append((edge, x_int, y_int))

        return intersections
