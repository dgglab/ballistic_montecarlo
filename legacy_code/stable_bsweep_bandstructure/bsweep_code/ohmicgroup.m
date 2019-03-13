classdef ohmicgroup < handle
    
    properties
        Nedges
        x_s
        y_s
        p_s
        L
        edgenorm
        edge_inject
        nedge_in_ohmic
    end
    
    methods
        function obj=ohmicgroup(frm_in,nedge_in_ohmic)
            obj.nedge_in_ohmic=nedge_in_ohmic;
            obj.Nedges=length(nedge_in_ohmic);
            obj.x_s=zeros(obj.Nedges,2);
            obj.y_s=zeros(obj.Nedges,2);
            obj.p_s=zeros(obj.Nedges,1);
            
            for i=1:obj.Nedges
               obj.x_s(i,1)=frm_in.px(nedge_in_ohmic(i)+1);
               obj.x_s(i,2)=frm_in.px(nedge_in_ohmic(i));
               obj.y_s(i,1)=frm_in.py(nedge_in_ohmic(i)+1);
               obj.y_s(i,2)=frm_in.py(nedge_in_ohmic(i));
            end
            
            L_s=sqrt((obj.x_s(:,2)-obj.x_s(:,1)).^2+(obj.y_s(:,2)-obj.y_s(:,1)).^2);
            obj.L=sum(L_s);
            obj.p_s=L_s/obj.L;
            
            obj.edgenorm=frm_in.edgenorm(nedge_in_ohmic);
            obj.edge_inject=0;
        end
        
        
        function new_edge_in(obj)
            p_temp=rand();
            p_cumulative=cumsum(obj.p_s);
            obj.edge_inject=find(p_temp<p_cumulative,1,'first');
             
        end
        
        function [x,y,inject_angle]=new_inject(obj)
           new_edge_in(obj);
                   
           inject_pos=rand();
           inject_angle=mod(randcos(-pi/2,pi/2)+obj.edgenorm(obj.edge_inject),2*pi);
    
            x=obj.x_s(obj.edge_inject,1)+inject_pos*(obj.x_s(obj.edge_inject,2)-obj.x_s(obj.edge_inject,1));
            y=obj.y_s(obj.edge_inject,1)+inject_pos*(obj.y_s(obj.edge_inject,2)-obj.y_s(obj.edge_inject,1));
            
        end
        
        function [x,y,norm]=inject_position(obj)
            new_edge_in(obj);
            
            inject_pos=rand();
            while inject_pos == 0 || inject_pos == 1
                % Don't inject exactly on end of edge
                inject_pos=rand();
            end
            
            x=obj.x_s(obj.edge_inject,1)+inject_pos*(obj.x_s(obj.edge_inject,2)-obj.x_s(obj.edge_inject,1));
            y=obj.y_s(obj.edge_inject,1)+inject_pos*(obj.y_s(obj.edge_inject,2)-obj.y_s(obj.edge_inject,1));
            norm = obj.edgenorm(obj.edge_inject);
            
        end
       
        function plot_random_inject(obj)
            [x,y,inject_angle]=new_inject(obj);
            
            ds=.05;
            dx=ds*cos(inject_angle);
            dy=ds*sin(inject_angle);
            
            plot(obj.x_s',obj.y_s','k',x,y,'r.',[x x+dx],[y y+dy],'r');
            
            
        end
    end
    
end

