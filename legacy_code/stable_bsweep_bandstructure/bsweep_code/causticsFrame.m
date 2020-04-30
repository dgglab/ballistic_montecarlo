classdef causticsFrame < handle
    properties
       px
       py
       edgestyle
       ms
       phis
       bs
       Ls
       edgenorm
       Nedges
       edgestats
       ogs 
       ohmics
       terminalOhmics
       injector_ohmic
    end
    methods
        function obj = causticsFrame(px,py,edgestyle)
            
            obj.px=px;
            obj.py=py;
            obj.edgestyle=edgestyle;
            
            obj.Nedges=length(px)-1;
            obj.ms=zeros(1,obj.Nedges);
            obj.phis=zeros(1,obj.Nedges);
            obj.bs=zeros(1,obj.Nedges);
            
            for i=1:obj.Nedges
                obj.ms(i)=(obj.py(i+1)-obj.py(i))/(obj.px(i+1)-obj.px(i));
                obj.phis(i)=atan2(obj.py(i+1)-obj.py(i),obj.px(i+1)-obj.px(i));
                obj.bs(i)=(obj.py(i)*obj.px(i+1)-obj.py(i+1)*obj.px(i))/(obj.px(i+1)-obj.px(i));
               % obj.edgenorm(i)=obj.phis(i)-pi/2+pi*inpolygon((obj.px(i+1)+obj.px(i))/2+obj.inPGtest*cos(obj.phis(i)+pi/2),(obj.py(i+1)+obj.py(i))/2+obj.inPGtest*sin(obj.phis(i)+pi/2),obj.px,obj.py);
                obj.edgenorm(i)=obj.phis(i)-pi/2;%+pi*inpolygon((obj.px(i+1)+obj.px(i))/2+obj.inPGtest*cos(obj.phis(i)+pi/2),(obj.py(i+1)+obj.py(i))/2+obj.inPGtest*sin(obj.phis(i)+pi/2),obj.px,obj.py);
                obj.Ls(i)=sqrt((obj.px(i+1)-obj.px(i)).^2+(obj.py(i+1)-obj.py(i)).^2);
            end
            
        end
        
        %% In Bounds?
        function out = inbounds(obj,x,y)
            out=inpolygon(x,y,obj.px,obj.py);
        end
        
        %% Line Crossed?
        function [xi, yi, crossedlines, corner] = linescrossed(obj,x1,x0,y1,y0)
           
            for n=1:obj.Nedges
                x13=obj.px(n)-x0;
                y13=obj.py(n)-y0;
                x21=obj.px(n+1)-obj.px(n);
                y21=obj.py(n+1)-obj.py(n);
                x43=x1-x0;
                y43=y1-y0;
                ts(n)=(x13*y21-x21*y13)/(x43*y21-x21*y43);
                us(n)=(x13*y43-x43*y13)/(x43*y21-x21*y43);
            end
            % Round to account for numerical error
            ts = round(ts,8);
            us = round(us,8);

            intersections=find(ts>=0 & ts<= 1 & ...
                us >= 0 & us <= 1);
            xi = [];
            yi = [];
            crossedlines = [];
            corner = [];

            
            
            if(length(intersections)>1)
                               
                [xi,yi,ii]=polyxpoly([x1 x0],[y1 y0],obj.px,obj.py);
                edgescrossed = ii(:,2);
                dLs=sqrt((x0-xi).^2+(y0-yi).^2);
                %[~, inds] = sort(dLs(dLs>0));
                [dLs, inds] = sort(dLs);
                %edgescrossed = edgescrossed(dLs>0);
                xi = xi(inds);
                yi = yi(inds);
                crossedlines=edgescrossed(inds);
                
                
                if length(intersections) ~= length(dLs)
                    % if missing a point, calculate distance and sort
                    extracross = setdiff(intersections,edgescrossed);
                    pxe = [obj.px(extracross+1) obj.px(extracross)];
                    pye = [obj.py(extracross+1) obj.py(extracross)];
                    [xic, yic] = polyxpoly([x1 x0], [y1 y0],pxe,pye);
                    
                    edgescrossed = [crossedlines; extracross'];
                    dLsc = sqrt((x0-xic).^2+(y0-yic).^2);
                    dLs = [dLs; dLsc];
                    xi = [xi; xic];
                    yi = [yi; yic];
                    [dLs, inds] = sort(dLs);
                    xi = xi(inds);
                    yi = yi(inds);
                    crossedlines=edgescrossed(inds);
                end
                
                % Detect if we hit a corner and store
                for i = 1:length(crossedlines)-1
                    if dLs(i) == dLs(i+1)
                        corner(end+1,:) = [crossedlines(i) crossedlines(i+1)];
                    end
                end
            elseif ~isempty(intersections)
                [xi,yi]=polyxpoly([x1 x0],[y1 y0],obj.px,obj.py);
                crossedlines = intersections;
            end
            
        end
        
        %% Specular Reflection
        function thetaout = specular_reflection(obj, thetain,crossedline)
            thetaout=2*obj.phis(crossedline)-thetain;
        end
        
        %% Nearly Specular Reflection
        % Simulates nearly specular reflection by adding a gaussian
        % distribution with speficied width (gausswidth) in theta.  Will
        % NOT ensure the reflected angle is still pointing out of the wall.
        function thetaout = nearly_specular_reflection(obj, thetain, crossedline, gausswidth)
            lineAngle = obj.phis(crossedline);
            thetaPerfect = 2*lineAngle - thetain;
            thetaout = thetaPerfect + gausswidth*randn(1,1);
        end
              
        %% Circle Cross
        function [xs, ys, linescrossed]=circle_cross(obj,c_arc)
            x2=obj.px(2:obj.Nedges+1)-c_arc.x0;
            x1=obj.px(1:obj.Nedges)-c_arc.x0;
            y2=obj.py(2:obj.Nedges+1)-c_arc.y0;
            y1=obj.py(1:obj.Nedges)-c_arc.y0;

            a=(x2-x1).^2+(y2-y1).^2;
            b=2*(x1.*x2-x1.*x1+y1.*y2-y1.*y1);
            c=x1.^2+y1.^2-c_arc.r^2;

            sqrtarg=sqrt(b.^2-4*a.*c);

            soln1=(-b+sqrtarg)./a/2;
            soln2=(-b-sqrtarg)./a/2;

            soln1_cross=imag(soln1)==0&real(soln1)<1&real(soln1)>0;
            soln2_cross=imag(soln2)==0&real(soln2)<1&real(soln2)>0;

            r1=find(soln1_cross);
            r2=find(soln2_cross);

            linescrossed=[r1,r2];
            xs=[x1(r1)+(x2(r1)-x1(r1)).*soln1(r1),x1(r2)+(x2(r2)-x1(r2)).*soln2(r2)]+c_arc.x0;
            ys=[y1(r1)+(y2(r1)-y1(r1)).*soln1(r1),y1(r2)+(y2(r2)-y1(r2)).*soln2(r2)]+c_arc.y0;

        end

            
        %% Nearest Intersection 
         function [xout,yout,edgeout]=nearest_intersection(obj,arc_in)
             [xs,ys,edge]=obj.circle_cross(arc_in);
            a=arc_in.arc_lengths(xs,ys);
            rmin=find(a==min(a));
            xout=xs(rmin);
            yout=ys(rmin);
            edgeout=edge(rmin);
            arc_in.L=a(rmin);
         end
        
        %% Next Nearest Intersection
        function [xout,yout,edgeout]=next_nearest_intersection(obj,arc_in)
            [xs,ys,edge]=obj.circle_cross(arc_in);
            a=arc_in.arc_lengths(xs,ys);
            mins=min(mod(a,abs(arc_in.r)*2*pi),mod(-a,abs(arc_in.r)*2*pi));
            a(find(mins==min(mins)))=inf;
            rmin=find(a==min(a));
            xout=xs(rmin);
            yout=ys(rmin);
            edgeout=edge(rmin);
            arc_in.L=a(rmin);
        end
        
        %% random reflection off the edge
        function theta_out=randCos_off_edge(obj,edge_in)
            theta_out=mod(obj.edgenorm(edge_in)+randcos(-pi/2,pi/2),2*pi);
        end
        
        %% Add Ohmics
        function obj = add_ohmics(obj,ohmics)
            
            x1=obj.px;
            y1=obj.py;
            obj.ohmics = ohmics;
            edgestyle1=zeros(1,size(obj.edgestyle,2));
            
            for i=1:length(ohmics)
                if ~isfield(obj.ohmics{i},'ohmNumber')
                    obj.ohmics{i}.ohmNumber=i;
                end
                [x1,y1,edgestyle1]=poly_overlap_cut(x1,y1,ohmics{i}.px,ohmics{i}.py,edgestyle1,obj.ohmics{i}.ohmNumber);
            end
            
            [x1,y1]=poly2cw(x1,y1);
            n_divider=find(isnan(x1));
            
            if(~isempty(n_divider))
                for i=1:length(n_divider)-1
                    [x1((n_divider(i-1)+1):(n_divider(i)-1)), y1((n_divider(i-1)+1):(n_divider(i)-1))]=poly2ccw(x1((n_divider(i-1)+1):(n_divider(i)-1)),y1((n_divider(i-1)+1):(n_divider(i)-1)));
                end
                [x1((n_divider(end)+1):length(x1)), y1((n_divider(end)+1):length(x1))]=poly2ccw(x1((n_divider(end)+1):length(x1)),y1((n_divider(end)+1):length(x1)));
            end
            obj=causticsFrame(x1,y1,edgestyle1);
            obj.ohmics = ohmics;
            obj.set_og();
        end
            
       
        %% Make Ohmics group
        function set_og(obj)
        clear obj.og
            for i=1:max(obj.edgestyle)
                    obj.ogs{i}=ohmicgroup(obj,find(obj.edgestyle==i));
           end 
        end
        

        %%
        function renumberOhmics(obj)
            maxNum=0;
            ohmicNumbers=zeros(1,length(obj.ohmics));
            for i=1:length(obj.ohmics)
                ohmicNumbers(i)=obj.ohmics{i}.ohmNumber;
            end
            
            i=1;
            while i<max(ohmicNumbers)
                if isempty(ohmicNumbers(ohmicNumbers==i))
                    ohmicNumbers(ohmicNumbers>i)=ohmicNumbers(ohmicNumbers>i)-1;
                else
                    i=i+1;
                end
            end
            for i=1:length(obj.ohmics)
                obj.ohmics{i}.ohmNumber=ohmicNumbers(i);
            end
        end
        
        %%
        function plotOhmicLabels(obj)
            plot(obj.px,obj.py);
            hold on
            for i=1:length(obj.ohmics)
                plot(obj.ohmics{i}.px,obj.ohmics{i}.py,'r');
                text(mean(obj.ohmics{i}.px(1:end-1)),mean(obj.ohmics{i}.py(1:end-1)),...
                    num2str(obj.ohmics{i}.ohmNumber),'FontWeight','Bold');
            end
            
            hold off
        end
        
        %%
        function out=terminal_ohmic_test(obj,edge)
            out=sum(obj.edgestyle(edge)==obj.terminalOhmics);
        end
        

        
        
        
        
        
        
        
        
    end
end