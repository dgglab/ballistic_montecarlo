import matplotlib.pyplot as plt
import ezdxf
from ezdxf.groupby import groupby
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from shapely.geometry import LineString


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
    edges: A list of tuples containing the line segement and edgestyle indicating the layer
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
                # ValueError is not technically the right error here
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

        xs, ys = self.body.exterior.coords.xy
        edgestyle = [0] * (len(xs) - 1)

        for layer in self.ps.keys():
            if layer == '0':
                continue
            for poly in self.ps[layer]:
                ohm = Polygon(poly)
                for i, (x, y) in enumerate(zip(xs[:-1], ys[:-1])):
                    if ohm.contains(LineString([(x, y), (xs[i+1], ys[i+1])])):
                        edgestyle[i] = int(layer)

        self.edges = []
        for i, (x, y) in enumerate(zip(xs[:-1], ys[:-1])):
            self.edges.append(([(x, y), (xs[i+1], ys[i+1])], edgestyle[i]))

    def gen_fig(self):
        '''
        Return a figure colored by the edge style
        '''
        fig = plt.figure()
        for edge in self.edges:
            x, y = list(zip(*edge[0]))
            plt.plot(x, y, color='C'+str(edge[1]))

        return fig
