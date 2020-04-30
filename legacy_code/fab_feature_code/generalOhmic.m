% generalOhmic is a class that generates a generic ohmic contact given a
% list of points. Additionally, a box with a specific pad is generated
% around the ohmic

classdef generalOhmic < fabFeatureGroup
    
    properties
       px
       py
       pad = 0.1

    end
    methods
         
        function obj = generalOhmic(px,py)
            
            
            obj.px=px;
            obj.py=py;
            
            generateOhmic(obj)
            
        end
        
        % Computes all of the coordinates based on the defined injector
        % properties.
        function generateOhmic(obj)
            
            % Create the ohmic
            inject1=fabFeature(obj.px,obj.py,0);
            obj.registerFeature(inject1,0,0,0);
            
            % Box around ohmic
            ohm1 = fabFeature(obj.px,obj.py,1);
            obj.registerFeature(ohm1,0,0,0);
            ohm1.grow(obj.pad);
        end
        
    end
    
end
            
            
            
            
            
            
            
            
            
            
           
