import matplotlib as mpl
import numpy as np

class Device:
    def __init__(self):
        self.features = []

    def register_feature(self, feature, delx, dely, theta):
        """
        Register a feature to the device. Allows one to shift the position of
        the feature in the coordinates of the device
        """
        feature.rotate_and_offset(delx,dely,theta)
        self.features.append(feature)


    def register_feature_group(self, feature_group, delx, dely, theta):
        for feature in feature_group.features:
            self.register_feature(feature, delx, dely, theta)


    def rotate_and_offset_device(self, delx, dely, theta):
        for feature in self.features:
            feature.rotate_and_offset(delx, dely, theta)


    def join_features(self):

    def plot():

    def DXF_output():

    def grow():



class Feature:
    def __init__(self, p, layer):
        self.update_shape(p)
        self.init_Tmatrix()
        self.set_layer(layer)


    def init_Tmatrix():
        self.Tmatrix = np.matlib.eye(3)


    def set_layer(self, layer):
        self.layer = layer


    def update_shape(self, p):
    """
    update coordinates to be clockwise
    """
        if p[0] == p[-1]:
            self.p_base = p
        else:
            p.append(p[0])
            self.p_base = p

        self.poly2cw()
        self.perform_transformation()


    def poly2cw(self):
        '''
        Checks the orienation of the polygon. If its verticies are oriented
        ccw, the following sum will be negative. The verticies are then
        reversed to be cw.
        '''
        orientation = 0

        for i, (px,py) in enumerate(self.p_base):
            px2 = self.p_base[i+1][0]
            py2 = self.p_base[i+1][1]
            orientation += (px2 - px) * (py2 + py)

        if orientation < 0:
            self.p_base = self.p_base[::-1]


    def perform_transformation(self):
        self.p = []
        for px, py in self.p_base:
            self.p.append(
                self.Tmatrix[0,0]*px + self.Tmatrix[0,1]*py + self.Tmatrix[0,2],
                self.Tmatrix[1,0]*px + self.Tmatrix[1,1]*py + self.Tmatrix[1,2]
            )

    def rotate_and_offset():

    def grow():
