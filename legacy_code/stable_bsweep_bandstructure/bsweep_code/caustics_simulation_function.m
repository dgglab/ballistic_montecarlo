function [edgenum,arcdata]=caustics_simulation_function_BS(frmgrp,rad_curv,N_inject,...
    p_ifbounce_then_scatter,p_bounce_off_ohmic,L_scatter)

frm1=frmgrp.frms{1};
injector_ohmic=frm1.injector_ohmic;


rad_curv0=rad_curv;


k=1;
%edgenum is array counting number of electrons interacting with an edge
edgenum=zeros(size(frm1.edgestyle,2),1);

n_in=1;
x1s=zeros(1E6,1);
y1s=zeros(1E6,1);
edges=zeros(1E6,1);
arc1=caustic_arc(1,1,1,1);
%Primary while loop
while(n_in<=N_inject)
    %new electron injected
    [x1,y1,theta1]=frm1.ogs{injector_ohmic}.new_inject();
    
    edge=frm1.ogs{injector_ohmic}.nedge_in_ohmic(frm1.ogs{injector_ohmic}.edge_inject);
    
    
    %store coordinates
    x1s(k)=x1;
    y1s(k)=y1;
    edges(k)=edge;
    
    %active_tradjectory is a flag that is made false when an electron is
    %absorbed by a grounded ohmic
    active_tradjectory=1;
    while active_tradjectory
        %define inital arc
        arc1(k)=caustic_arc(x1,y1,rad_curv,theta1);
        
        %edge is zero if scattered in bulk. Below if-then statement finds
        %the first edge the caustic will hit. It distinguishes between
        %arcs started from edges and from the bulk, to prevent miscounting
        %the source edge.
        
        % [x1a,y1a,edgea,frmouta]=frmgrp.next_nearest_intersection(arc1(k));
        [x1a,y1a,edgea,frmouta]=frmgrp.nearest_intersection(arc1(k));
        if(abs(x1-x1a)<1E-6 & abs(y1-y1a)<1E-6)
            [x1a,y1a,edgea,frmouta]=frmgrp.next_nearest_intersection(arc1(k));
        end
        x1=x1a;y1=y1a;edge=edgea;frmout=frmouta;
        
        
        if(~isempty(x1))
            %if there's an edge, set a to the arc length between points
            a=arc1(k).arc_lengths(x1,y1);
            
            x1s(k)=x1;
            y1s(k)=y1;
            edges(k)=edge;
            
        else
            %arc length is infinite if no edge
            a=inf;
        end
        
        %L_nonstop is an exponential random variable to implement bulk
        %scattering
        L_nonstop=exprnd(L_scatter,1,1);
        
        
        if(L_nonstop<a)
            %scattering occurs if L_nonstop<a. sets a new x1,y1, and
            %theta1 for the next arc.
            [x1,y1,theta1]=arc1(k).move_L(L_nonstop);
            theta1=mod(theta1+randcos(-pi,pi),2*pi);
            
            x1s(k)=x1;
            y1s(k)=y1;
            arc1(k).setL(L_nonstop);
            edge=0;
            edges(k)=edge;
            k=k+1;
            
        else
            %if an edge is hit, check what class of frame it comes from
            switch class(frmout)
                case 'causticsFrame'
                    
                    if ~frmout.edgestyle(edge)
                        %case of non-ohmic edge
                        if rand<p_ifbounce_then_scatter
                            theta1=frmout.randCos_off_edge(edge);
                        else
                            theta1=frmout.specular_reflection(arc1(k).incident_angle(x1,y1),edge);
                        end
                        %theta1=frmout.randCos_off_edge(edge);
                        k=k+1;
                    elseif frmout.terminal_ohmic_test(edge)
                        
                        
                        if rand<p_bounce_off_ohmic
                            if rand<p_ifbounce_then_scatter
                                theta1=frmout.randCos_off_edge(edge);
                            else
                                theta1=frmout.specular_reflection(arc1(k).incident_angle(x1,y1),edge);
                            end
                            k=k+1;
                        else
                            %sunk to ground
                            active_tradjectory=0;
                        end
                        
                        
                    else
                        
                        if rand<p_bounce_off_ohmic
                            if rand<p_ifbounce_then_scatter
                                theta1=frmout.randCos_off_edge(edge);
                            else
                                theta1=frmout.specular_reflection(arc1(k).incident_angle(x1,y1),edge);
                            end
                            k=k+1;
                        else
                            %absorbs and reinjects
                            edgenum(edge)=edgenum(edge)+1;
                            [x1,y1,theta1]=frmout.ogs{frm1.edgestyle(edge)}.new_inject();
                            k=k+1;
                        end
                    end
                    
                case 'causticsBarrier'
                    %refracts or reflects off of barrier
                    [theta1,rad_curv]=frmout.reflect_or_transmit(arc1(k).incident_angle(x1,y1),rad_curv0,edge);
                    k=k+1;
                    
            end
        end
        
    end
    k=k+1;
    
    edgenum(edge)=edgenum(edge)+1;
    if(~sum(edge==frm1.ogs{injector_ohmic}.nedge_in_ohmic))
        n_in=n_in+1;
    end
end


Narcs=length(arc1);

arcdata=zeros(5,Narcs);
for i=1:Narcs
    arcdata(:,i)=[arc1(i).x1,arc1(i).y1,-arc1(i).r,arc1(i).theta0,arc1(i).L];
end

