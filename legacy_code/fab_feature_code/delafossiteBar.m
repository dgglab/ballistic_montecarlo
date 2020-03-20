% collimatingHallbar is a class that creats a Hall bar with four
% collimating injectors.

classdef delafossiteBar < fabFeatureGroup
    
    properties
        deviceWidth
        deviceLength
        
    end
    
    methods
        
        function obj = delafossiteBar(lengthBar,widthBar)
            
            obj.deviceLength = lengthBar;
            obj.deviceWidth = widthBar;
            
            generateDevice(obj)
            
        end
        
        function generateDevice(obj)
            
            % Define the device body
            W = obj.deviceWidth;
            L = obj.deviceLength;
                        
            px=[-W/2 W/2 W/2 -W/2];
            py=[L/2 L/2 -L/2 -L/2];
            
            body0=fabFeature(px,py,0);
            obj.registerFeature(body0,0,0,0);
            
        end

        
        
        % Adds collecting ohmics on the sides
        function addOhmics(obj)
            W = obj.deviceWidth;
            L = obj.deviceLength;
            pad = 0.1;
            
            
            % Bottom Ohmic
            px = [-W/2 - pad, W/2 + pad, W/2 + pad, -W/2 - pad];
            py = [-L/2 + pad, -L/2 + pad, -L/2 - pad, -L/2 - pad];
            ohm0 = fabFeature(px,py,1);
            obj.registerFeature(ohm0,0,0,0);
                       
            % Top Ohmic
            px = [-W/2 - pad, W/2 + pad, W/2 + pad, -W/2 - pad];
            py = [L/2 + pad, L/2 + pad, L/2 - pad, L/2 - pad];
            ohm1 = fabFeature(px,py,1);
            obj.registerFeature(ohm1,0,0,0);
            
        end

    end
    
end











