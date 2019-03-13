classdef caustic_arc < handle
    
    properties
       x1
       y1
       x0
       y0
       r
       omega0
       theta0
       from_line=0;
       L=0;
    end
    
    methods
       
        function obj=caustic_arc(x1,y1,r,theta)
            obj.x1=x1;
            obj.y1=y1;
            obj.r =-r;
            obj.theta0=theta;
            obj.x0=obj.x1-obj.r*sin(theta);
            obj.y0=obj.y1+obj.r*cos(theta);
            obj.omega0=mod(atan2(obj.y1-obj.y0,obj.x1-obj.x0),2*pi);
        
        end
        
        
        function omega=angle_meas(obj,x2,y2)
            x2 = real(x2);
            y2 = real(y2);
            x10=obj.x1-obj.x0;
            y10=obj.y1-obj.y0;
            x20=x2-obj.x0;
            y20=y2-obj.y0;
            x12=obj.x1-x2;
            y12=obj.y1-y2;
            
            L20=sqrt(x20^2+y20^2);
            L12=sqrt(x12^2+y12^2);
            
            omega=asin((x10*y20-x20*y10)/obj.r/L20);
            
            if L12^2>(obj.r^2+L20^2)
                omega=sign(omega)*(pi-abs(omega));
            end
            omega=mod(omega,2*pi);
            %end
            %omega = 0;
            
        end
        
        function theta1=incident_angle(obj,x1,y1)
            theta1=mod(obj.angle_meas(x1,y1)*sign(obj.r)+obj.theta0,2*pi);
        end
        
        function [xout,yout,thetaout]=move_L(obj,Lin)
            omega=Lin/obj.r+obj.omega0;
            xout=abs(obj.r)*cos(omega)+obj.x0;
            yout=abs(obj.r)*sin(omega)+obj.y0;
            thetaout=mod(obj.theta0+Lin/obj.r,2*pi);
            
        end
        
        function plot_arc(obj,x2,y2)
            omega=angle_meas(obj,x2,y2);
            
            omegas=linspace(0,sign(obj.r)*omega,200);
            
            %plot([obj.x0,obj.x1,x2],[obj.y0,obj.y1,y2],'ro',obj.x0+abs(obj.r)*cos(omegas+obj.omega0),obj.y0+abs(obj.r)*sin(omegas+obj.omega0));
            plot([obj.x1,x2],[obj.y1,y2],'ro',obj.x0+abs(obj.r)*cos(omegas+obj.omega0),obj.y0+abs(obj.r)*sin(omegas+obj.omega0));
            
            
        end
        
        function [xout,yout]=arc_points(obj,L,dL)
            N=round(L/dL);
            omega=L/abs(obj.r);
            omegas=linspace(0,sign(obj.r)*omega,N);
            xout=obj.x0+abs(obj.r)*cos(omegas+obj.omega0);
            yout=obj.y0+abs(obj.r)*sin(omegas+obj.omega0);
            
            
        end
        
        function L_out=arc_lengths(obj,xs,ys)
           N=length(xs);
           L_out=zeros(N,1);
           for i=1:N
              L_out(i)=abs(obj.r)*angle_meas(obj,xs(i),ys(i)); 
           end
           
        end
        
        %afunction [xout,yout]=nearest_intersection(obj,xs,ys)
        %    a=obj.arc_lengths(xs,ys);
        %    xout=xs(find(a==min(a)));
        %    yout=ys(find(a==min(a)));
        %end
        
        
        
         function [xout,yout,edgeout]=nearest_intersection(obj,frm)
             [xs,ys,edge]=frm.circle_cross(obj);
            a=obj.arc_lengths(xs,ys);
            rmin=find(a==min(a));
            xout=xs(rmin);
            yout=ys(rmin);
            edgeout=edge(rmin);
            obj.L=a(rmin);
        end
        
        function [xout,yout,edgeout]=next_nearest_intersection(obj,frm)
            [xs,ys,edge]=frm.circle_cross(obj);
            a=obj.arc_lengths(xs,ys);
            mins=min(mod(a,abs(obj.r)*2*pi),mod(-a,abs(obj.r)*2*pi));
            a(find(mins==min(mins)))=inf;
            rmin=find(a==min(a));
            xout=xs(rmin);
            yout=ys(rmin);
            edgeout=edge(rmin);
            obj.L=a(rmin);
        end
        function setL(obj,Lin)
            obj.L=Lin;
        end
        
    end
    
    
end