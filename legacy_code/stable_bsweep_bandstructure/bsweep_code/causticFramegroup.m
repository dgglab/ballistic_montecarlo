classdef  causticFramegroup < handle
   
    properties
        frms
        N_frms
        Nedges_in_frms
        edge_offsets
        px
        py
        super_frm
        edgestyles
    end
    
    
    methods
        function obj=causticFramegroup
            obj.N_frms=0;
            obj.px=[];
            obj.py=[];
            obj.Nedges_in_frms=[];
        end
        
        
        
        function register_frame(obj,frame_in)
            obj.N_frms=obj.N_frms+1;
            obj.frms{obj.N_frms}=frame_in;
            obj.Nedges_in_frms=[obj.Nedges_in_frms,length(frame_in.px)-1];
            if obj.N_frms==1
                obj.edge_offsets=0;
                obj.edgestyles=obj.N_frms*ones(1,length(frame_in.px)-1);
                obj.px=[frame_in.px];
                obj.py=[frame_in.py];
            else
                
                obj.px=[obj.px, nan, frame_in.px];
                obj.py=[obj.py, nan, frame_in.py];
                obj.edge_offsets=[obj.edge_offsets,obj.edge_offsets(end)+obj.Nedges_in_frms(end-1)+2];
                obj.edgestyles=[obj.edgestyles,nan,nan,obj.N_frms*ones(1,length(frame_in.px)-1)];
            end
            obj.super_frm=causticsFrame(obj.px,obj.py,obj.edgestyles);
        end
        
        function overwrite_last_frame(obj, frame_in)
            obj.frms{obj.N_frms} = frame_in;
            
            obj.Nedges_in_frms = [obj.Nedges_in_frms(1:end-1), length(frame_in.px)-1];
            
            last_frame_index = find(isnan(obj.px), 1, 'last');
            obj.px = [obj.px(1:last_frame_index), frame_in.px];
            obj.py = [obj.py(1:last_frame_index), frame_in.py];
            last_frame_index = find(isnan(obj.edgestyles), 1, 'last');
            obj.edgestyles=[obj.edgestyles(1:last_frame_index),obj.N_frms*ones(1,length(frame_in.px)-1)];
            obj.edge_offsets = [obj.edge_offsets(1:end-1), obj.edge_offsets(end)+obj.Nedges_in_frms(end-1)+2];
            
        end
        
        
        
        function plotFrmLabels(obj)
            for i=1:length(obj.frms)
                if isa(obj.frms{i},'causticsFrame')
                    plot(obj.frms{i}.px,obj.frms{i}.py);
                else
                    plot(obj.frms{i}.px,obj.frms{i}.py,'r');
                end
                hold on;
                text(mean(obj.frms{i}.px(1:end-1)),mean(obj.frms{i}.py(1:end-1)),...
                    num2str(i),'FontWeight','Bold');
            end
            
            hold off;
        end
                
        
        function setCausticsBarrier(obj,frmNum,E_rel)
           if isa(obj.frms{frmNum},'causticsBarrier')
               
               obj.frms{frmNum}=causticsBarrier(obj.frms{frmNum}.px,...
                   obj.frms{frmNum}.py,E_rel);
           else
               sprintf('not a causticsBarrier!\n');
           end
        end
        
        function    [xout,yout,edgeout,frmout]=next_nearest_intersection(obj,arc_in)

            [xout,yout,edgeout1]=obj.super_frm.next_nearest_intersection(arc_in);
            frmnum=obj.super_frm.edgestyle(edgeout1);
            edgeout=edgeout1-obj.edge_offsets(frmnum);
            frmout=obj.frms{frmnum};

        end
        
        
        function    [xout,yout,edgeout,frmout]=nearest_intersection(obj,arc_in)

            [xout,yout,edgeout1]=obj.super_frm.nearest_intersection(arc_in);
            if edgeout1>0
                frmnum=obj.super_frm.edgestyle(edgeout1);
                edgeout=edgeout1-obj.edge_offsets(frmnum);
                frmout=obj.frms{frmnum};
            else
                edgeout=0;
                frmout=obj.frms{1};
            end

        end
        
        function [xi, yi, edgeout,corner,frmout] = linescrossed(obj,x1,x0,y1,y0)
            
            [xi, yi, edgeout1, corner] = obj.super_frm.linescrossed(x1,x0,y1,y0);
            if edgeout1>0
                frmnum=obj.super_frm.edgestyle(edgeout1);
                
                for i = 1:size(frmnum,2)
                    edgeout(i)=edgeout1(i)-obj.edge_offsets(frmnum(i));
                    frmout{i}=obj.frms{frmnum(i)};
                end
            else
                edgeout=0;
                frmout=obj.frms{1};
            end
           
%             for n=1:obj.Nedges_in_frms
%                 x13=obj.px(n)-x0;
%                 y13=obj.py(n)-y0;
%                 x21=obj.px(n+1)-obj.px(n);
%                 y21=obj.py(n+1)-obj.py(n);
%                 x43=x1-x0;
%                 y43=y1-y0;
%                 ts(n)=(x13*y21-x21*y13)/(x43*y21-x21*y43);
%                 us(n)=(x13*y43-x43*y13)/(x43*y21-x21*y43);
%             end
%             
%             crossedline=find(ts>0&ts<1&us>0&us<1);
%             if(length(crossedline)>1)
%                 [xi,yi,ii]=polyxpoly([x1 x0],[y1 y0],obj.px,obj.py);
%                 dLs=sqrt((x0-xi).^2+(y0-yi).^2);
%                 crossedline=ii(dLs==min(dLs));
%             end
            
        end
        
        function out = cornerCheck(obj,corner)
            % cornerCheck takes the cross product of the edges forming
            % a corner to determine if the corner is 'interior/exterior'
            
            edge = max(corner);
               
            P0 = [obj.px(edge-1), obj.py(edge-1)];
            P1 = [obj.px(edge), obj.py(edge)];
            P2 = [obj.px(edge+1), obj.py(edge+1)];
            
            % Normalized vectors
            n1 = (P1 - P0) / norm(P1 - P0);  
            n2 = (P2 - P1) / norm(P2 - P1);
            
            out = det([n2;n1]);
            
        end

    end
          
        
    
    
    
    
    
    
end