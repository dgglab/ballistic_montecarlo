% collimatingHallbar is a class that creats a Hall bar with four
% collimating injectors.

classdef delafossite_design_2 < fabFeatureGroup
    
    properties
        dr1Width
        dr1Pos
        dr2Width
        dr2Pos
        inPos
        deviceLength
        deviceWidth
        
    end
    
    methods
        
        function obj = delafossite_design_2(deviceLength,deviceWidth)
            
            obj.deviceLength = deviceLength;
            obj.deviceWidth = deviceWidth;
            
            generateDevice(obj)
            
        end
        
        function generateDevice(obj)
            
            % Define the device body
            W = obj.deviceWidth;
            L = obj.deviceLength;
                        
            px=[-L/2 L/2 L/2 -L/2];
            py=[W/2 W/2 -W/2 -W/2];
            
            body1=fabFeature(px,py,0);
            obj.registerFeature(body1,0,0,0);
            
        end
        
        % Adds collimators
        function addDrains(obj,dr1Width,dr1Pos,...
                dr2Width,dr2Pos)
            
            obj.dr1Width = dr1Width;
            obj.dr1Pos = dr1Pos;
            obj.dr2Width = dr2Width;
            obj.dr2Pos = dr2Pos;
            d = 5;
                        
            % Define the drains
            % Left drain
            W = obj.deviceWidth;
            L = obj.deviceLength;
            px =  [-L/2 + obj.dr1Pos, -L/2 + obj.dr1Pos + obj.dr1Width,...
                -L/2 + obj.dr1Pos + obj.dr1Width, -L/2 + obj.dr1Pos];
            py = [W/2, W/2, W/2 + d, W/2 + d];
            
            dr1 = fabFeature(px,py,0);
            obj.registerFeature(dr1,0,0,0);
            
            % Right drain
            px =  [-L/2 + obj.dr2Pos, -L/2 + obj.dr2Pos + obj.dr2Width,...
                -L/2 + obj.dr2Pos + obj.dr2Width, -L/2 + obj.dr2Pos];
            py = [W/2, W/2, W/2 + d, W/2 + d];
            
            dr2 = fabFeature(px,py,0);
            obj.registerFeature(dr2,0,0,0);
            
            
            % Add the ohmics
            % Left ohmic
            pad = 0.1;
            px =  [-L/2 + obj.dr1Pos - pad, -L/2 + obj.dr1Pos + obj.dr1Width + pad,...
                -L/2 + obj.dr1Pos + obj.dr1Width + pad, -L/2 + obj.dr1Pos - pad];
            py = [W/2 + d - pad, W/2 + d - pad, W/2 + d + pad, W/2 + d + pad];
            
            ohm1 = fabFeature(px,py,1);
            obj.registerFeature(ohm1,0,0,0);
            
            % Right ohmic
            px =  [-L/2 + obj.dr2Pos - pad, -L/2 + obj.dr2Pos + obj.dr2Width + pad,...
                -L/2 + obj.dr2Pos + obj.dr2Width + pad, -L/2 + obj.dr2Pos - pad];
            py = [W/2 + d - pad, W/2 + d - pad, W/2 + d + pad, W/2 + d + pad];
            
            ohm2 = fabFeature(px,py,1);
            obj.registerFeature(ohm2,0,0,0);
            
            
        end
        
        
        % Adds collecting ohmics on the sides
        function addInjectors(obj,inPos)
            W = obj.deviceWidth;
            L = obj.deviceLength;
            obj.inPos = -L/2 + inPos;
            inApp = 0.25;
            inCone = inApp*1.5;
            coneL = 3 * inApp;
            inSep = 1;
            %traceL = 3.75;
            %slope = 1/6;
            pad = 0.1;
            
            % in0
            n = 0;
            px = n*1 + [obj.inPos, obj.inPos - inCone/2, obj.inPos + inApp + inCone/2, obj.inPos + inApp];
            py = [-W/2, -W/2 - coneL, -W/2 - coneL, -W/2];
                                 
            in0 = fabFeature(px,py,0);
            obj.registerFeature(in0,0,0,0);
            
%             % Trace 0
%             px = n*1 + [obj.inPos + inApp + inCone/2, obj.inPos - inCone/2, ...
%                 obj.inPos - inCone/2 - traceL*slope, obj.inPos + inApp + inCone/2 - traceL*slope];
%             
%             py = [-W/2 - coneL, -W/2 - coneL, -W/2 - coneL - traceL, -W/2 - coneL- traceL];
%             
%             tr0 = fabFeature(px,py,0);
%             obj.registerFeature(tr0,0,0,0);
            
            % Ohmic 0
            px = n*1 + [obj.inPos - inCone/2 - pad, obj.inPos - inCone/2 - pad, ...
                obj.inPos + inApp + inCone/2 + pad, obj.inPos + inApp + inCone/2 + pad];
            py = [-W/2 - coneL + pad, -W/2 - coneL - pad,...
                -W/2 - coneL - pad, -W/2 - coneL + pad];
            ohmIn0 = fabFeature(px,py,1);
            obj.registerFeature(ohmIn0,0,0,0);
            
            
            % in1
            n = 1;
            px = n*1 + [obj.inPos, obj.inPos - inCone/2, obj.inPos + inApp + inCone/2, obj.inPos + inApp];
            py = [-W/2, -W/2 - coneL, -W/2 - coneL, -W/2];
                                 
            in1 = fabFeature(px,py,0);
            obj.registerFeature(in1,0,0,0);
            
%             % Trace 1
%             px = n*1 + [obj.inPos + inApp + inCone/2, obj.inPos - inCone/2, ...
%                 obj.inPos - inCone/2 - traceL*slope, obj.inPos + inApp + inCone/2 - traceL*slope];
%             
%             py = [-W/2 - coneL, -W/2 - coneL, -W/2 - coneL - traceL, -W/2 - coneL- traceL];
%             
%             tr1 = fabFeature(px,py,0);
%             obj.registerFeature(tr1,0,0,0);
            
            % Ohmic 1
            px = n*1 + [obj.inPos - inCone/2 - pad, obj.inPos - inCone/2 - pad, ...
                obj.inPos + inApp + inCone/2 + pad, obj.inPos + inApp + inCone/2 + pad];
            py = [-W/2 - coneL + pad, -W/2 - coneL - pad,...
                -W/2 - coneL - pad, -W/2 - coneL + pad];
            ohmIn1 = fabFeature(px,py,1);
            obj.registerFeature(ohmIn1,0,0,0);
            
            % in2
            n = 2;
            px = n*1 + [obj.inPos, obj.inPos - inCone/2, obj.inPos + inApp + inCone/2, obj.inPos + inApp];
            py = [-W/2, -W/2 - coneL, -W/2 - coneL, -W/2];
                                 
            in2 = fabFeature(px,py,0);
            obj.registerFeature(in2,0,0,0);
            
%             % Trace 2
%             px = n*1 + [obj.inPos + inApp + inCone/2, obj.inPos - inCone/2, ...
%                 obj.inPos - inCone/2 - traceL*slope, obj.inPos + inApp + inCone/2];
%             
%             py = [-W/2 - coneL, -W/2 - coneL, -W/2 - coneL - traceL, -W/2 - coneL- traceL];
%             
%             tr2 = fabFeature(px,py,0);
%             obj.registerFeature(tr2,0,0,0);
            
            % Ohmic 2
            px = n*1 + [obj.inPos - inCone/2 - pad, obj.inPos - inCone/2 - pad, ...
                obj.inPos + inApp + inCone/2 + pad, obj.inPos + inApp + inCone/2 + pad];
            py = [-W/2 - coneL + pad, -W/2 - coneL - pad,...
                -W/2 - coneL - pad, -W/2 - coneL + pad];
            ohmIn2 = fabFeature(px,py,1);
            obj.registerFeature(ohmIn2,0,0,0);
            
            
            % in3
            n = 3;
            px = n*1 + [obj.inPos, obj.inPos - inCone/2, obj.inPos + inApp + inCone/2, obj.inPos + inApp];
            py = [-W/2, -W/2 - coneL, -W/2 - coneL, -W/2];
                                 
            in3 = fabFeature(px,py,0);
            obj.registerFeature(in3,0,0,0);
            
            % Ohmic 3
            px = n*1 + [obj.inPos - inCone/2 - pad, obj.inPos - inCone/2 - pad, ...
                obj.inPos + inApp + inCone/2 + pad, obj.inPos + inApp + inCone/2 + pad];
            py = [-W/2 - coneL + pad, -W/2 - coneL - pad,...
                -W/2 - coneL - pad, -W/2 - coneL + pad];
            ohmIn3 = fabFeature(px,py,1);
            obj.registerFeature(ohmIn3,0,0,0);
            
%             % Trace 3
%             
%             px = n*1 + [obj.inPos + inApp + inCone/2, obj.inPos - inCone/2, ...
%                 obj.inPos - inCone/2, obj.inPos + inApp + inCone/2 + traceL*slope];
%             
%             py = [-W/2 - coneL, -W/2 - coneL, -W/2 - coneL - traceL, -W/2 - coneL- traceL];
%             
%             tr3 = fabFeature(px,py,0);
%             obj.registerFeature(tr3,0,0,0);
            
            
            % in4
            n = 4;
            px = n*1 + [obj.inPos, obj.inPos - inCone/2, obj.inPos + inApp + inCone/2, obj.inPos + inApp];
            py = [-W/2, -W/2 - coneL, -W/2 - coneL, -W/2];
                                 
            in4 = fabFeature(px,py,0);
            obj.registerFeature(in4,0,0,0);
            
%             % Trace 4
%             
%             px = n*1 + [obj.inPos + inApp + inCone/2, obj.inPos - inCone/2, ...
%                 obj.inPos - inCone/2 + traceL*slope, obj.inPos + inApp + inCone/2 + traceL*slope];
%             
%             py = [-W/2 - coneL, -W/2 - coneL, -W/2 - coneL - traceL, -W/2 - coneL- traceL];
%             
%             tr4 = fabFeature(px,py,0);
%             obj.registerFeature(tr4,0,0,0);
            
            % Ohmic 4
            px = n*1 + [obj.inPos - inCone/2 - pad, obj.inPos - inCone/2 - pad, ...
                obj.inPos + inApp + inCone/2 + pad, obj.inPos + inApp + inCone/2 + pad];
            py = [-W/2 - coneL + pad, -W/2 - coneL - pad,...
                -W/2 - coneL - pad, -W/2 - coneL + pad];
            ohmIn4 = fabFeature(px,py,1);
            obj.registerFeature(ohmIn4,0,0,0);
            
            % in5
            n = 5;
            px = n*1 + [obj.inPos, obj.inPos - inCone/2, obj.inPos + inApp + inCone/2, obj.inPos + inApp];
            py = [-W/2, -W/2 - coneL, -W/2 - coneL, -W/2];
                                 
            in5 = fabFeature(px,py,0);
            obj.registerFeature(in5,0,0,0);
            
%             % Trace 5
%             
%             px = n*1 + [obj.inPos + inApp + inCone/2, obj.inPos - inCone/2, ...
%                 obj.inPos - inCone/2 + traceL*slope, obj.inPos + inApp + inCone/2 + traceL*slope];
%             
%             py = [-W/2 - coneL, -W/2 - coneL, -W/2 - coneL - traceL, -W/2 - coneL- traceL];
%             
%             tr5 = fabFeature(px,py,0);
%             obj.registerFeature(tr5,0,0,0);
            
            % Ohmic 5
            px = n*1 + [obj.inPos - inCone/2 - pad, obj.inPos - inCone/2 - pad, ...
                obj.inPos + inApp + inCone/2 + pad, obj.inPos + inApp + inCone/2 + pad];
            py = [-W/2 - coneL + pad, -W/2 - coneL - pad,...
                -W/2 - coneL - pad, -W/2 - coneL + pad];
            ohmIn5 = fabFeature(px,py,1);
            obj.registerFeature(ohmIn5,0,0,0);
            
            
        end

    end
    
end











