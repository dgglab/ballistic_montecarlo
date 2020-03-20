% collimatingHallbar is a class that creats a Hall bar with four
% collimating injectors.

classdef hallBar < fabFeatureGroup
    
    properties
        deviceWidth
        probeWidth
        probeLength
        bodySQ
        injSQ
        
    end
    
    methods
        
        function obj = hallBar(width,bodySQ,injSQ)
            
            obj.deviceWidth = width;
            obj.probeWidth = 0.3;
            obj.probeLength = 0.5;
            obj.bodySQ = bodySQ;
            obj.injSQ = injSQ;

            
            generateTopGate(obj)
            
        end
        
        function generateTopGate(obj)
            
            % Define the device body
            W = obj.deviceWidth;
            ohmW = obj.probeWidth;
            ohmL = obj.probeLength;
            bodyL = W*obj.bodySQ;
            injL = W*obj.injSQ;
                        
            px=[-bodyL/2-ohmW-injL -bodyL/2-ohmW -bodyL/2-ohmW -bodyL/2 -bodyL/2 ...
                bodyL/2 bodyL/2 bodyL/2+ohmW bodyL/2+ohmW bodyL/2+ohmW+injL ...
                bodyL/2+ohmW+injL bodyL/2+ohmW bodyL/2+ohmW bodyL/2 bodyL/2 ...
                -bodyL/2 -bodyL/2 -bodyL/2-ohmW -bodyL/2-ohmW -bodyL/2-ohmW-injL ...
                -bodyL/2-ohmW-injL];
            py=[W/2 W/2 W/2+ohmL W/2+ohmL W/2 ...
               W/2 W/2+ohmL W/2+ohmL W/2 W/2 ...
               -W/2 -W/2 -W/2-ohmL -W/2-ohmL -W/2 ...
               -W/2 -W/2-ohmL -W/2-ohmL -W/2 -W/2 ...
               W/2];
            
            body0=fabFeature(px,py,0);
            obj.registerFeature(body0,0,0,0);
            
        end

        
        
        % Adds collecting ohmics on the sides
        function addEtch(obj)
            W = obj.deviceWidth;
            ohmW = obj.probeWidth;
            ohmL = obj.probeLength;
            bodyL = W*obj.bodySQ;
            injL = W*obj.injSQ;
            pad = 0.1;
            
            
            px=[-bodyL/2-ohmW-injL+pad -bodyL/2-ohmW+pad -bodyL/2-ohmW+pad -bodyL/2-pad -bodyL/2-pad ...
                bodyL/2+pad bodyL/2+pad bodyL/2+ohmW-pad bodyL/2+ohmW-pad bodyL/2+ohmW+injL-pad ...
                bodyL/2+ohmW+injL-pad bodyL/2+ohmW-pad bodyL/2+ohmW-pad bodyL/2+pad bodyL/2+pad ...
                -bodyL/2-pad -bodyL/2-pad -bodyL/2-ohmW+pad -bodyL/2-ohmW+pad -bodyL/2-ohmW-injL+pad ...
                -bodyL/2-ohmW-injL+pad];
            py=[W/2-pad W/2-pad W/2+ohmL-pad W/2+ohmL-pad W/2-pad ...
               W/2-pad W/2+ohmL-pad W/2+ohmL-pad W/2-pad W/2-pad ...
               -W/2+pad -W/2+pad -W/2-ohmL+pad -W/2-ohmL+pad -W/2+pad ...
               -W/2+pad -W/2-ohmL+pad -W/2-ohmL+pad -W/2+pad -W/2+pad ...
               W/2-pad];
            ohm0 = fabFeature(px,py,1);
            obj.registerFeature(ohm0,0,0,0);
                       
            
        end

    end
    
end











