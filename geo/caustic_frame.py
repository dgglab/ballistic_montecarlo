import numpy as np
import matplotlib.pyplot as plt
import ezdxf
from ezdxf.groupby import groupby
from shapely.geometry import Polygon
from shapely.geometry import LineString
from shapely.geometry.polygon import orient


class Frame:
    '''
    Frame to used in ballistic monte carlo
    Required layers:
        0: Body of device, must be a single polygon
        1: Injecting contact(s)
        2: Grounded contacts(s)

    Optional layers:
        3-n: floating contacts, grouped by layer

    Attributes
    ps: dictionary containing the verticies of all polygons in the device separated by layer
    body: shapely polygon of the body segmented according to all ohmics
    edges: A list of tuples containing the line segement and layer indicating the layer
           of the edges
    '''

    def __init__(self, dxf_path):
        '''
        loads in a dxf file and generates a frame from the file
        '''
        doc = ezdxf.readfile(dxf_path)
        layers = groupby(entities=doc.modelspace(), dxfattrib='layer')

        if len(layers['0']) != 1:
            if len(layers['0']) == 0:
                raise ValueError(
                    'Body (layer 0) does not exist, assign to layer 0')
            else:
                raise ValueError('Body (layer 0) is disjoint')

        if len(layers['1']) == 0:
            raise ValueError('No injector (layer 1)')

        if len(layers['2']) == 0:
            raise ValueError('No ground (layer 2)')

        self._extract_points(layers)
        self._gen_frame()
        self._gen_matrices_for_det()

    def _extract_points(self, layers):
        '''
        extracts all the points of polygons into a dictionary
        Input:
            layers: a grouped imported dxf file
        '''
        self.ps = {}
        for layer in layers.keys():
            for poly in layers[layer]:
                poly_ps = []
                for vert in list(poly):
                    poly_ps.append((vert[0], vert[1]))

                if layer in self.ps:
                    self.ps[layer].append(poly_ps)
                else:
                    self.ps[layer] = [poly_ps]

    def _gen_frame(self):
        '''
        generate the frame by intersecting the contacts with the body
        '''
        for layer in self.ps.keys():
            if layer == '0':
                self.body = Polygon(self.ps[layer][0])
            else:
                for poly in self.ps[layer]:
                    ohm = Polygon(poly)
                    inter = self.body.intersection(ohm)
                    self.body = self.body.union(inter)

        self.body = orient(self.body, -1)  # order clockwise
        xs, ys = self.body.exterior.coords.xy
        layers = [0] * (len(xs) - 1)

        for layer in self.ps.keys():
            if layer == '0':
                continue
            for poly in self.ps[layer]:
                ohm = Polygon(poly)
                for i, (x, y) in enumerate(zip(xs[:-1], ys[:-1])):
                    if ohm.contains(LineString([(x, y), (xs[i+1], ys[i+1])])):
                        layers[i] = int(layer)

        self.edges = []
        for i, (x, y) in enumerate(zip(xs[:-1], ys[:-1])):
            self.edges.append(Edge(x, y, xs[i+1], ys[i+1], layers[i]))

    def _gen_matrices_for_det(self):
        '''
        Generates a few matricies that will be used in calculating intersections with edges
        '''
        self.px0 = np.array([edge.xs[0] for edge in self.edges])
        self.px1 = np.array([edge.xs[1] for edge in self.edges])
        self.py0 = np.array([edge.ys[0] for edge in self.edges])
        self.py1 = np.array([edge.ys[1] for edge in self.edges])
        self.x23 = self.px0 - self.px1
        self.y23 = self.py0 - self.py1

    def gen_fig(self):
        '''
        Return a figure colored by the edge style
        '''
        fig = plt.figure()
        for edge in self.edges:
            x, y = [edge.xs, edge.ys]
            if edge.layer == 0:
                plt.plot(x, y, color='k')
            else:
                plt.plot(x, y, color='C'+str(edge.layer-1))

        return fig

    def get_inject_position(self, layer, bias=True):
        '''
        Randomly pick a position along the perimeter of a given ohmic edge
        '''
        edges_in_layer = []
        for edge in self.edges:
            if edge.layer == layer:
                edges_in_layer.append(edge)

        if not edges_in_layer:
            raise ValueError('No edges in given layer')

        r = 0  # We are forcing r in (0, 1), to not return a corner
        while r == 0:
            r = np.random.rand()

        lengths = [edge.length for edge in edges_in_layer]
        edge_prob = np.cumsum(lengths)/np.sum(lengths)
        ind = np.argmax(edge_prob > r)

        injecting_edge = edges_in_layer[ind]

        # Scale the random number to the fractional position along the chosen edge
        if ind > 0:
            r_scaled = (r - edge_prob[ind - 1])/(lengths[ind]/np.sum(lengths))
        else:
            r_scaled = r/(lengths[ind]/np.sum(lengths))

        (x_inject, y_inject) = injecting_edge.interpolate_frac_positon(r_scaled)

        if bias:
            bias_vector = 1E-10 * \
                np.array([np.cos(injecting_edge.normal_angle),
                          np.sin(injecting_edge.normal_angle)])
            x_inject += bias_vector[0]
            y_inject += bias_vector[1]

        return (x_inject, y_inject), injecting_edge


class Edge:
    '''
    Edge represents a given edge in a frame. It stores properties of the edge, such
    as its coordiantes and the probability distribtution for injection from that edge
    '''

    def __init__(self, x0, y0, x1, y1, layer):
        self.start = (x0, y0)
        self.end = (x1, y1)
        self.xs = [x0, x1]
        self.ys = [y0, y1]
        self.layer = layer
        self.length = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        self.normal_angle = np.arctan2(y1-y0, x1-x0) - np.pi/2
        self.normal = np.array(
            [np.cos(self.normal_angle), np.sin(self.normal_angle)])
        self.linestring = LineString([(x0, y0), (x1, y1)])

    def __repr__(self):
        return '(({0}, {1}), ({2}, {3}), {4})'.format(self.start[0], self.start[1], self.end[0], self.end[1], self.layer)

    def interpolate_frac_positon(self, f):
        '''
        Linarly interpolate a position on the edge according to a given fractional
        '''
        x = self.start[0] + f * (self.end[0] - self.start[0])
        y = self.start[1] + f * (self.end[1] - self.start[1])
        return(x, y)

    def set_in_prob(self, in_prob, cum_prob):
        self.in_prob = in_prob
        self.cum_prob = cum_prob

    def get_injection_index(self):
        return Edge.compute_injection_index(self.cum_prob)

    @staticmethod
    def compute_injection_index(cum_prob):
        '''
        Returns an index n which is allowable to inject into and a fractional scaling factor
        '''
        r = np.random.rand()
        n = np.argmax(cum_prob > r)
        if n == 0:
            f = 1 - (r) / (cum_prob[n])
        else:
            f = 1 - (r - cum_prob[n-1])/(cum_prob[n] - cum_prob[n-1])
        return (n, f)
