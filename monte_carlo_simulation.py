import time
import pickle
import os.path

import numpy as np
from shapely.geometry import Point
from enum import IntEnum
from copy import deepcopy

from ballistic_montecarlo.geo.caustic_frame import Edge, OhmicLines
from ballistic_montecarlo.bandstructure.caustic_bandstructure import Bandstructure


def calc_ohmstats(fields, results):
    ohmstats = {}
    zeros = np.zeros(np.shape(fields))
    for i, (r, _) in enumerate(zip(results, fields)):
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


ALL_STATES = set(range(11))


class Simulation:
    def __init__(self, frame, k, phi, field, p_scatter=1.0, p_ohmic_absorb=1.0, ohmic_lines=OhmicLines([])):
        '''
        inputs:
            frame: instance of a Frame representing the device to simulate
            k: list of arrays wave vectors defining Fermi surface [kx, ky]
            phi: angle of crystal axis relative to the device
            field: magnetic field to simulate
            p_scatter: probability of scattering from an edge (instead of reflecting)
            p_ohmic_absorb: probability of an ohmic absorbing an impinging charge
        '''
        self._frame = deepcopy(frame)
        self._ohmic_lines = deepcopy(ohmic_lines)
        self._bandstructure = Bandstructure(k, phi, field)
        self._p_scatter = p_scatter
        self._p_ohmic_absorb = p_ohmic_absorb

        # Calculate injection probabilities for each edge.
        max_layer = 0
        for edge in self._frame.edges:
            in_prob, cum_prob = self._bandstructure.calculate_injection_prob(
                edge.normal_angle)
            edge.set_in_prob(in_prob, cum_prob)
            if edge.layer > max_layer:
                max_layer = edge.layer

        for line in self._ohmic_lines.lines_as_edges:
            line.layer += max_layer + 1

    def run_simulation_with_cache(self, identifier, n_inject, path='data/', stored_states=ALL_STATES, debug=False):
        full_path = path+identifier+'.pkl'
        if os.path.isfile(path):
            print('path {} already exists, loading data'.format(full_path))
            with open(full_path, 'rb') as f:
                edge_to_count, trajectories = pickle.load(f)
            return edge_to_count, trajectories

        print('path {} does not exist, running simulation'.format(full_path))
        edge_to_count, trajectories = self.run_simulation(
            n_inject, stored_states=stored_states, debug=False)
        with open(full_path, 'wb') as f:
            pickle.dump((edge_to_count, trajectories), f)

        return edge_to_count, trajectories

    def run_simulation(self, n_inject, stored_states=ALL_STATES, debug=False):
        '''
        Propagates n_inject charge carriers until they are absorbed by a grounded contact
        '''

        trajectories = []
        edge_to_count = {}
        for edge in self._frame.edges:
            edge_to_count[edge] = 0

        for line in self._ohmic_lines.lines_as_edges:
            edge_to_count[line] = 0

        for _ in range(n_inject):
            trajectory = []
            state = TrajectoryState.INJECTING
            (x, y), injecting_edge = self._frame.get_inject_position(1)
            n_f = injecting_edge.get_injection_index()

            step_params = [(n_f, x, y, state, injecting_edge)]

            # Optimize to either or not based on full or empty filter
            trajectory.extend(
                filter(lambda step: step[-2] in stored_states, step_params))

            while state != TrajectoryState.ABSORBED:
                n_f = step_params[-1][0]
                x = step_params[-1][1]
                y = step_params[-1][2]
                step_params, line_crosses = self._step_position(
                    n_f, x, y, debug)

                state = step_params[-1][-2]

                trajectory.extend(
                    filter(lambda step: step[-2] in stored_states, step_params))

                for step_param in step_params:
                    edge = step_param[-1]
                    if edge in edge_to_count:
                        edge_to_count[edge] += 1

                for line in line_crosses:
                    edge_to_count[line] += 1

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
        line_crosses = self._get_crosses(step_coords)

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
                    n_f_new = self._specular(n_f_int, edge)

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

            if layer == 1:
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

        return step_params, line_crosses

    def _update_position(self, n_f_in, x_in, y_in):
        [x_out, y_out] = [x_in, y_in] + n_f_in[1] * \
            self._bandstructure.dr[:, n_f_in[0]]
        n_f_out = ((n_f_in[0] + 1) % (np.shape(self._bandstructure.dr)[1]), 1)
        return n_f_out, x_out, y_out

    def _specular(self, n_f, edge):
        xd = edge.xs[0] - edge.xs[1] 
        yd = edge.ys[0] - edge.ys[1] 
        #get the slope of the edge
        line_slope = float('inf') 
        if xd != 0:
            line_slope = yd/xd
        line_x1 = self._bandstructure.r[0][n_f[0]]
        line_y1 = self._bandstructure.r[1][n_f[0]]
        matching_segments = list()
        # go through all the segments and try to find the intersecting point
        #   if an intersection exists, add it to matching_segments
        for real_segment in range(0, len(self._bandstructure.r[0])):
            next_real_segment = (real_segment + 1) % len(self._bandstructure.r[0])
            segment_x3 = self._bandstructure.r[0][real_segment]
            segment_y3 = self._bandstructure.r[1][real_segment]
            segment_x4 = self._bandstructure.r[0][next_real_segment]
            segment_y4 = self._bandstructure.r[1][next_real_segment]
            if line_slope == float('inf'):
                if (line_x1 >= segment_x3 and line_x1 < segment_x4) or (line_x1 > segment_x4 and line_x1 <= segment_x3):
                    matching_segments.append(real_segment)
            else:
                line_x2 = 100.0 #random number to get another point
                line_y2 = ((line_x2 - line_x1) * line_slope) / line_y1
                denominator = (line_x1 - line_x2)*(segment_y3 - segment_y4) - (line_y1 - line_y2)*(segment_x3 - segment_x4)
                if denominator !=0:
                    intersection_x = ((line_x1*line_y2 - line_y1*line_x2) * (segment_x3 - segment_x4) - (line_x1 - line_x2) * (segment_x3*segment_y4 - segment_y3*segment_x4)) / denominator
                    intersection_y = ((line_x1*line_y2 - line_y1*line_x2) * (segment_y3 - segment_y4) - (line_y1 - line_y2) * (segment_x3*segment_y4 - segment_y3*segment_x4)) / denominator
                    if (intersection_x <= segment_x3 and intersection_x >= segment_x4) or (intersection_x <= segment_x4 and intersection_x >= segment_x3):
                        matching_segments.append(real_segment)
        #we find the one that matches the current position and then we discard every other one as
        # those would get out of the figure
        new_matching_segments = list()
        for existing_finder in range(len(matching_segments)):
            if matching_segments[existing_finder] == n_f[0]:
                goes_out_of_the_figure = True
                for odd_detector in range(len(matching_segments)):
                    if not goes_out_of_the_figure:
                        new_matching_segments.append(matching_segments[(odd_detector+existing_finder)%len(matching_segments)])
                    goes_out_of_the_figure = not goes_out_of_the_figure 
        #pythagoras on each segment to check which one is the closest
        min_point = -1
        min_distance = float('inf')
        for testing_point in new_matching_segments:
            segment_x3 = self._bandstructure.r[0][testing_point]
            segment_y3 = self._bandstructure.r[1][testing_point]
            this_distance = np.sqrt(np.power(line_x1-segment_x3, 2)+(np.power(line_y1-segment_y3, 2)))
            if min_distance > this_distance:
                min_distance = this_distance
                min_distance = testing_point
        if min_distance == -1:
            die('Can\'t find an intersection for specular reflection')
            return
        else:
            return (min_point, 1)

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
        ts, us = self.get_ts_us(x01, x02, self._frame.x23,
                                y01, y02, self._frame.y23)

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

    def get_ts_us(self, x01, x02, x23, y01, y02, y23):
        with np.errstate(divide='ignore'):
            ts = (x02*y23 - y02*x23) / (x01*y23 - y01*x23)
            us = -(x01*y02 - y01*x02) / (x01*y23 - y01*x23)
        return ts, us

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

    def _get_crosses(self, step_coords):
        x = step_coords[0][0]
        y = step_coords[0][1]
        x_new = step_coords[1][0]
        y_new = step_coords[1][1]

        x_del = x_new - x
        y_del = y_new - y
        x01 = -x_del
        y01 = -y_del
        x02 = x - self._ohmic_lines.px0
        y02 = y - self._ohmic_lines.py0

        crosses = []
        ts, us = self.get_ts_us(x01, x02, self._ohmic_lines.x23,
                                y01, y02, self._ohmic_lines.y23)
        for i, (t, u) in enumerate(zip(ts, us)):
            if 0 <= t and t <= 1 and 0 <= u and u <= 1:
                line = self._ohmic_lines.lines_as_edges[i]
                crosses.append(line)
        return crosses
