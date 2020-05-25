import time

import numpy as np
from shapely.geometry import Point
from enum import IntEnum

from geo.caustic_frame import Edge
from bandstructure.caustic_bandstructure import Bandstructure


def calc_ohmstats(fields, results):
    ohmstats = {}
    zeros = np.zeros(np.shape(fields))
    for i, (r, field) in enumerate(zip(results, fields)):
        for edge in r.get()[0].keys():
            count = r.get()[0][edge]
            if edge.layer in ohmstats:
                ohmstats[edge.layer][i] += count
            else:
                ohmstats[edge.layer] = zeros.copy()
                ohmstats[edge.layer][i] += count
    return ohmstats


class TrajectoryState(IntEnum):
    '''
    Used for specifying the state of a step returned by Simulation._step_positionS
    '''
    INJECTING = 1
    PROPAGATE = 2
    COLLISION = 3
    SCATTER = 4
    REFLECT = 5
    ABSORBED = 6
    CCOLLISION = 7
    CSCATTER = 8
    CREFLECT = 9
    CABSORBED = 10
    ERROR = 11


class Simulation:
    def __init__(self, frame, k, phi, field, p_scatter=1.0, p_ohmic_absorb=1.0):
        '''
        inputs:
            frame: instance of a Frame representing the device to simulate
            k: list of arrays wave vectors defining Fermi surface [kx, ky]
            phi: angle of crystal axis relative to the device
            field: magnetic field to simulate
            p_scatter: probability of scattering from an edge (instead of reflecting)
            p_ohmic_absorb: probability of an ohmic absorbing an impinging charge
        '''
        self._frame = frame
        self._bandstructure = Bandstructure(k, phi, field)
        self._p_scatter = p_scatter
        self._p_ohmic_absorb = p_ohmic_absorb

        # Calculate injection probabilities for each edge.
        for edge in self._frame.edges:
            in_prob, cum_prob = self._bandstructure.calculate_injection_prob(
                edge.normal_angle)
            edge.set_in_prob(in_prob, cum_prob)

    def run_simulation(self, n_inject, debug=False):
        '''
        Propagates n_inject charge carriers until they are absorbed by a grounded contact
        '''
        trajectories = []
        edge_to_count = {}
        for edge in self._frame.edges:
            edge_to_count[edge] = 0

        for _ in range(n_inject):
            trajectory = []
            state = TrajectoryState.INJECTING
            (x, y), injecting_edge = self._frame.get_inject_position(1)
            n_f = injecting_edge.get_injection_index()

            step_params = [(n_f, x, y, state, injecting_edge)]
            trajectory.append(step_params[0])

            while state != TrajectoryState.ABSORBED:
                n_f = step_params[-1][0]
                x = step_params[-1][1]
                y = step_params[-1][2]
                step_params = self._step_position(n_f, x, y, debug)

                state = step_params[-1][-2]
                trajectory.extend(step_params)
                for step_param in step_params:
                    edge = step_param[-1]
                    if edge in edge_to_count:
                        edge_to_count[edge] += 1

            trajectories.append(trajectory)

        return edge_to_count, trajectories

    def _step_position(self, n_f, x, y, debug=False):
        step_params = []
        if debug and not self._frame.body.intersects(Point(x, y)):
            print('Previous step stepped out of bounds')
            print(n_f, x, y)
            return [(n_f, x, y, TrajectoryState.ERROR, None)]

        n_f_new, x_new, y_new = self._update_position(n_f, x, y)

        step_coords = ([(x, y), (x_new, y_new)])

        intersections = self._get_sorted_intersections(step_coords)

        # TODO: This big ole block is convoluted and has lots of repetition.
        #       Consider refactoring, or at least naming things more clearly.
        if len(intersections) == 0:
            step_params.append(
                (n_f_new, x_new, y_new, TrajectoryState.PROPAGATE, None))

        elif len(intersections) == 1 or intersections[0][3] != intersections[1][3]:
            # Single edge intersection
            edge, x_int, y_int, _ = intersections[0]
            n_f_int = self._get_n_f_intersection(
                n_f, [(x, y), (x_int, y_int), (x_new, y_new)])

            if edge.layer == 0:
                # Device edge
                step_params.append(
                    (n_f_int, x_int, y_int, TrajectoryState.COLLISION, edge))

                if np.random.rand() < self._p_scatter:
                    n_f_new = self._scatter(edge)
                    step_params.append(
                        (n_f_new, x_int, y_int, TrajectoryState.SCATTER, None))
                else:
                    # Replace with specular reflection
                    n_f_new = self._scatter(edge)
                    step_params.append(
                        (n_f_new, x_int, y_int, TrajectoryState.REFLECT, None))

            elif edge.layer == 2:
                # Grounded ohmic
                if np.random.rand() < self._p_ohmic_absorb:
                    # Absorbed
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.ABSORBED, edge))
                else:
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.COLLISION, None))

                    if np.random.rand() < self._p_scatter:
                        n_f_new = self._scatter(edge)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.SCATTER, None))
                    else:
                        # Replace with specular reflection
                        n_f_new = self._scatter(edge)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.REFLECT, None))
            else:
                # Generic ohmic
                if np.random.rand() < self._p_ohmic_absorb:
                    # Absorb and reemit
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.ABSORBED, edge))

                    (x_new, y_new), reinjecting_edge = self._frame.get_inject_position(
                        edge.layer)
                    n_f_new = reinjecting_edge.get_injection_index()

                    step_params.append(
                        (n_f_new, x_new, y_new, TrajectoryState.INJECTING, None))
                else:
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.COLLISION, None))

                    if np.random.rand() < self._p_scatter:
                        n_f_new = self._scatter(edge)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.SCATTER, None))
                    else:
                        # Replace with specular reflection
                        n_f_new = self._scatter(edge)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.REFLECT, None))
        else:
            # Corner intersection
            edge_0, x_int, y_int, _ = intersections[0]
            edge_1, _, _, _ = intersections[1]

            n_f_int = self._get_n_f_intersection(
                n_f, [(x, y), (x_int, y_int), (x_new, y_new)])

            edges = [edge_0, edge_1]
            index = np.argmax([edge_0.layer, edge_1.layer])
            layer = edges[index].layer
            edge_for_count = edges[index]  # to avoid double counting

            if layer == 0:
                # Device edge
                step_params.append(
                    (n_f_int, x_int, y_int, TrajectoryState.CCOLLISION, edge_for_count))

                if np.random.rand() < self._p_scatter:
                    n_f_new = self._corner_scatter(edge_0, edge_1)
                    step_params.append(
                        (n_f_new, x_int, y_int, TrajectoryState.CSCATTER, None))
                else:
                    # Replace with specular reflection
                    n_f_new = self._corner_scatter(edge_0, edge_1)
                    step_params.append(
                        (n_f_new, x_int, y_int, TrajectoryState.CREFLECT, None))

            elif layer == 2:
                # Grounded ohmic
                if np.random.rand() < self._p_ohmic_absorb:
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.ABSORBED, edge_for_count))
                else:
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.CCOLLISION, None))

                    if np.random.rand() < self._p_scatter:
                        n_f_new = self._corner_scatter(edge_0, edge_1)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.CSCATTER, None))
                    else:
                        # Replace with specular reflection
                        n_f_new = self._corner_scatter(edge_0, edge_1)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.CREFLECT, None))
            else:
                # Generic ohmic
                if np.random.rand() < self._p_ohmic_absorb:
                    # Absorb and reemit
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.CABSORBED, edge_for_count))

                    (x_new, y_new), reinjecting_edge = self._frame.get_inject_position(
                        layer)
                    n_f_new = reinjecting_edge.get_injection_index()

                    step_params.append(
                        (n_f_new, x_new, y_new, TrajectoryState.INJECTING, None))

                else:
                    step_params.append(
                        (n_f_int, x_int, y_int, TrajectoryState.CCOLLISION, None))

                    if np.random.rand() < self._p_scatter:
                        n_f_new = self._corner_scatter(edge_0, edge_1)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.CSCATTER, None))
                    else:
                        # Replace with specular reflection
                        n_f_new = self._corner_scatter(edge_0, edge_1)
                        step_params.append(
                            (n_f_new, x_int, y_int, TrajectoryState.CREFLECT, None))

        return step_params

    def _update_position(self, n_f_in, x_in, y_in):
        [x_out, y_out] = [x_in, y_in] + n_f_in[1] * \
            self._bandstructure.dr[:, n_f_in[0]]
        n_f_out = ((n_f_in[0] + 1) % (np.shape(self._bandstructure.dr)[1]), 1)
        return n_f_out, x_out, y_out

    # Is this method really necessary?
    def _scatter(self, edge):
        return edge.get_injection_index()

    def _corner_scatter(self, edge_0, edge_1):
        # Convolve the probility distributions
        in_prob = (edge_0.in_prob*edge_1.in_prob) / \
            np.sum(edge_0.in_prob*edge_1.in_prob)
        cum_prob = np.cumsum(in_prob)
        return Edge.compute_injection_index(cum_prob)

    def _get_sorted_intersections(self, step_coords):
        x = step_coords[0][0]
        y = step_coords[0][1]
        intersections = self._get_intersections(step_coords)

        def compute_s(intersection):
            edge, x_int, y_int = intersection
            S = np.sqrt((x_int - x)**2 + (y_int - y)**2)
            return (edge, x_int, y_int, S)

        intersections_with_S = map(compute_s, intersections)
        if len(intersections) < 2:
            return list(intersections_with_S)
        else:
            return sorted(intersections_with_S, key=lambda intersections: intersections[3])

    def _get_intersections(self, step_coords, bias=True):
        '''
        Given a step, check all the edges in the _frame to see if the step crosses
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
        x02 = x - self._frame.px0
        y02 = y - self._frame.py0

        intersections = []

        ts = (x02*self._frame.y23 - y02*self._frame.x23) / \
             (x01*self._frame.y23 - y01*self._frame.x23)
        us = -(x01*y02 - y01*x02) / (x01*self._frame.y23 - y01*self._frame.x23)

        for i, (t, u) in enumerate(zip(ts, us)):
            if 0 <= t and t <= 1 and 0 <= u and u <= 1:
                edge = self._frame.edges[i]
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

    def _get_n_f_intersection(self, n_f, coords):
        '''
        Returns the scaling factor required to step from (x, y) to (x_int, y_int) given n_f
        '''
        x = coords[0][0]
        y = coords[0][1]
        x_int = coords[1][0]
        y_int = coords[1][1]
        x_new = coords[2][0]
        y_new = coords[2][1]
        f = np.sqrt(((x_int - x)**2 + (y_int-y)**2) /
                    ((x_new - x)**2 + (y_new - y)**2))
        return (n_f[0], f)
