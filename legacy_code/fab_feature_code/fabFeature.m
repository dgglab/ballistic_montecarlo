% fabFeature is a class structure for programattically defining fabrication
% features. It handles rotations and translations and establishes the layer
% number of the given polygon

classdef fabFeature < handle
    
    properties
        px_base,py_base,px,py
        Tmatrix
        N_edges
        LayerNum
    end
    
    methods
        
        %object constructor. 
        function obj = fabFeature(px,py,LayerNum)
            
            obj.initTmatrix();
            obj.updateShape(px,py);
            obj.setLayerNum(LayerNum);
            
        end
        
        %initialized Tmatrix for constructor and resetting coordinates
        function initTmatrix(obj)
            obj.Tmatrix=eye(3);
        end
        
        %sets the LayerNum of the fabFeature
        function setLayerNum(obj,LayerNum)
            obj.LayerNum=LayerNum;
        end
        
        %updates the base coordinates, and preserves the 
        function updateShape(obj,px,py)
            [px,py]=poly2cw(px,py);
            if (py(1)==py(end))&&(px(1)==px(end))||(sum(isnan(py))>0)
                obj.px_base=px;
                obj.py_base=py;
            else
                obj.px_base=[px,px(1)];
                obj.py_base=[py,py(1)];
            end
            obj.N_edges=length(obj.px_base);
            obj.performTransform();
            
        end
        
        %updates Tmatrix and performs the transform for px,py
        function rotate_and_offset(obj,theta,delx,dely)
            
            %Matrix Multiplication of previous transform and new transform
            obj.Tmatrix=[cos(theta),-sin(theta),delx;...
                     sin(theta), cos(theta),dely;...
                     0         , 0         , 1]*obj.Tmatrix;
            %calculate px,py from base and Tmatrix
            obj.performTransform();

        end
        
        %caculates px and py from base units and transform matrix
        function performTransform(obj)
            obj.px=obj.Tmatrix(1,1)*obj.px_base+...
                   obj.Tmatrix(1,2)*obj.py_base+...
                   obj.Tmatrix(1,3);
            obj.py=obj.Tmatrix(2,1)*obj.px_base+...
                   obj.Tmatrix(2,2)*obj.py_base+...
                   obj.Tmatrix(2,3);
        end
        
        %simple plot function
        function plot(obj,c_in)
            plot(obj.px,obj.py,c_in);
            axis equal;
        end
           
        % Grows or shrinks a given polygon
        function grow(obj,pad)
                        
            % figure(); hold on
            % plot(obj.px_base,obj.py_base)
            for i=1:size(obj.px_base,2)-1
                phis(i)=atan2(obj.py_base(i+1)-obj.py_base(i),obj.px_base(i+1)-obj.px_base(i));   
                edgenorm(i)=phis(i)+pi/2;
                
            end
            
            for i=1:size(edgenorm,2)
                ip = mod(i,size(obj.px_base,2))+1;
                dx = pad*cos(edgenorm(i));
                dy = pad*sin(edgenorm(i));
                edgex(i,:) = [obj.px_base(i) + dx obj.px_base(ip)+dx];
                edgey(i,:) = [obj.py_base(i) + dy obj.py_base(ip)+dy];
                
            end
            
            for i=1:size(edgenorm,2)
                ip = mod(i,size(edgenorm,2))+1;
                x1 = edgex(i,1); x2 = edgex(i,2);
                y1 = edgey(i,1); y2 = edgey(i,2);
                x3 = edgex(ip,1); x4 = edgex(ip,2);
                y3 = edgey(ip,1); y4 = edgey(ip,2);
                
                obj.px_base(ip) = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
                obj.py_base(ip) = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
            end
            
            obj.px_base(end) = obj.px_base(1);
            obj.py_base(end) = obj.py_base(1);
            obj.updateShape(obj.px_base,obj.py_base)
            %plot(obj.px_base,obj.py_base)
        end
            

    end
end
