from device_geo import Device,Feature

class Bar(Device):
    def __init__(self, bar_length, bar_width):
        self.bar_length = bar_length
        self.bar_width = bar_width

        self.generate_body()
        self.add_ohmics()
        self.join_features()
        self.plot()


    def generate_body(self):
        L = self.bar_length
        W = self.bar_width

        p = [(-W/2.0,  L/2.0),
             ( W/2.0,  L/2.0),
             ( W/2.0, -L/2.0),
             (-W/2.0, -L/2.0)]

        body0 = Feature(p, 0);
        self.register_feature(body0, 0, 0, 0);


    def add_ohmics(self):
        L = self.bar_length
        W = self.bar_width
        pad = 0.1

        # Bottom ohmic
        p = [(-W/2.0 - pad, -L/2.0 + pad),
             ( W/2.0 + pad, -L/2.0 + pad),
             ( W/2.0 + pad, -L/2.0 - pad),
             (-W/2.0 - pad, -L/2.0 - pad)]

        ohm0 = Feature(p, 1)
        self.register_feature(ohm0, 0, 0, 0)

        # Top ohmic
        p = [(-W/2.0 - pad, L/2.0 + pad),
             ( W/2.0 + pad, L/2.0 + pad),
             ( W/2.0 + pad, L/2.0 - pad),
             (-W/2.0 - pad, L/2.0 - pad)]

        ohm1 = Feature(p, 1)
        self.register_feature(ohm1, 0, 0, 0)
