function [BS,edgenum,arcdata]=caustics_simulation_function_BS(frmgrp,BSp,BSn,B_cur,N_inject,...
    p_ifbounce_then_scatter,p_bounce_off_ohmic,L_scatter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% caustics_simulation_function_BS.m
% Aaron Sharpe
% 10/17
%
% This program propigates electrons through time steps
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug = false; % Toggle debug plotting

frm1=frmgrp.frms{1};
injector_ohmic=frm1.injector_ohmic;


B_cur0=B_cur;

% Find the real space orbit for the given field
if B_cur0>0
    BS = BSp;
else
    BS = BSn;
end
BS.KtoR(B_cur0); % Need to change for parallel


k=1;
% Edgenum is array counting number of electrons interacting with an edge
% used in computing voltages at contacts
edgenum = zeros(size(frm1.edgestyle,2),1);

n_in = 1;             % Counter for injected electrons
x1s = zeros(1,1E6);   % Store fin/init points of arcs
y1s = zeros(1,1E6);   %
ind1s = zeros(2,1E6); % Store indices
edges = zeros(1,1E6); % Store edges seen



%Primary while loop
while(n_in<=N_inject)
    % New electron injected
    [x1,y1,norm] = frm1.ogs{injector_ohmic}.inject_position();
    ind1 = BS.inject(norm);
    

    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debuging
    %                                     x1 = 5.9999;
    %                                     y1 = 9.3411;
    %                                     ind1 = 62;
    %                             lastedge = 3;
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Find which edge the electron was injected from
    edge=frm1.ogs{injector_ohmic}.nedge_in_ohmic(frm1.ogs{injector_ohmic}.edge_inject);
    
    
    % Store initial coordinates
    x1s(k) = x1;
    y1s(k) = y1;
    x2 = x1; y2 = y1;
    ind2 = ind1;
    ind1s(1,k) = ind1;
    edges(k) = edge;
    
    kstore = k; % store k at start of current electron
    k=k+1; % Increment k
    
    
    % Store edge for decerning new edge hit
    lastedge = edge;
    lastEdgeReset = 0;
    lastcol = 0; % Number o steps since last collision
    
    
    % Active_tradjectory is a flag that is made false when an electron is
    % absorbed by a grounded ohmic
    active_trajectory = 1;
    
    error = 0; % Error flag for using try/catch
    
    edge = 0; % Initialize edge check
    
    justreflect = 1; % Just reflect flag
    
    % Storage of parameters
    store = [x1 y1 ind1 edge];
    sstore = store;
    ssstore = sstore;
    
    
    
    while active_trajectory
         try
            ssstore = sstore;
            sstore = store;
            store = [x1 y1 ind1 edge];
            
            
            inds = ind1;
            xs = x1; ys = y1;
            
            x1 = x2; y1 = y2; % Update coordinates
            ind1 = mod(ind2-1,size(BS.Rx,2)-1)+1; % Update index
            ind2 = floor(mod(ind1,size(BS.Rx,2)-1)+1); % Step to next integer
            
            % Compute new coordinates
            x2 = x1 + BS.Rx(ind2) - BS.fracInd(BS.Rx,ind1);
            y2 = y1 + BS.Ry(ind2) - BS.fracInd(BS.Ry,ind1);
            
            x2 = round(x2,4);
            y2 = round(y2,4);
            
            % Check inbounds here
%             if justreflect
%                 if ~frm1.inbounds(x2,y2) && ~frm1.inbounds(x1,y1)
%                     fprintf('inbounds check \n')
%                     figure(1); hold on
%                     plot(store(1),store(2),'.')
%                     plot(sstore(1),sstore(2),'.')
%                     plot(ssstore(1),ssstore(2),'.')
%                     drawnow
%                     ssstore
%                     sstore
%                     store
%                     
% %                     figure; hold on;
% %                     plot(BS.Rx,BS.Ry,'.-')
% %                     plot(BS.Rx(ind1r),BS.Ry(ind1r),'ok')
% %                     plot(BS.Rx(floor(ind2r)),BS.Ry(floor(ind2r)),'ob')
% %                     plot(BS.Rx(ceil(ind2r)),BS.Ry(ceil(ind2r)),'or')
% %                     edgenorm = frmout.edgenorm(edge);
% %                     testedge = find(edgenorm == BS.unorms);
% %                     BS.plotInProb(testedge,10)
%                     debug = true;
%                 end
%             end
            
            justreflect = 0;
            
            if debug
                figure(1); hold on
                plot(x1,y1,'.')
                plot(x2,y2,'.')
                drawnow
            end
            
            % If the number of steps since the last collision is larger
            % than the number of steps in the Fermi surface
            if lastcol > size(BS.Rx,2)
                if debug == false
                    fprintf('electron is stuck\n')
                end
                display([ind1 ind2 x1 x2 y1 y2])
%                 debug = true;
                fpritf('force error\n') % force an error
            end
            
            
            % Find all edges crossed sorted by distance from (x1,y1)
            [xi,yi,edge,corner,frmsout] = frmgrp.linescrossed(x2,x1,y2,y1);
            
            % Take the first edge different from the edge
            % that was crossed in the last time step
            
            if edge ~= 0 % If we saw an edge
                % Check for new edges crossed
                if ~isempty(corner)
                    % Corner handling
                    [crow, ~] = find(corner == lastedge);
                    
                    if all(ismember(corner(crow,:),lastedge))
                        % If we just hit a corner
                        
                        newedges = setdiff(edge,corner(crow,:),'stable');
                        corner(crow,:) = [];
                        
                    elseif any(ismember(corner(crow,:),lastedge))
                        % If we got rounded to a corner
                        out = frmgrp.cornerCheck(corner(crow,:));
                        if out >= 0 % Exterior corner or flat
                            newedges = setdiff(edge,lastedge,'stable');
                        else % Interior corner
                            newedges = setdiff(edge,corner(crow,:),'stable');
                        end
                        lastedge = corner(crow,:);
                        corner(crow,:) = [];
                        lastEdgeReset = 1;
                    end
                    
                else
                    % Just find new edges
                    newedges = setdiff(edge,lastedge,'stable');
                end
                
                if isempty(newedges)
                    % If no new nedges seen
                    edge = 0;
                else
                    if ~isempty(xi) % Able to calculate intersection
                        eind = find(edge==newedges(1),1);
                        edge = edge(eind);
                        xi = round(xi(eind),4);
                        yi = round(yi(eind),4);
                        
                        % find corner edge as well or clear corner
                        [crow, ~] = find(corner == edge);
                        if ~isempty(crow)
                            corner = corner(crow,:);
                        else
                            corner = [];
                        end
                        
                    else % Didnt actually detect intersection
                        edge = 0;
                        %debug = true;
                    end
                end
            end
            
            if isempty(edge)
                edge = 0;
            end
            

            if edge == 0
                % If no edge is seen
                
                lastedge = 0;
                lastEdgeReset = 0;
                lastcol = lastcol + 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                %
                %     Insert bulk scattering code
                %
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            else
                                         
                frmout = frmsout{eind};
				% Interpolate the index of the intersection
                [indf] = BS.intersectInd(x1,xi,x2,y1,yi,y2,ind1);
                indi = floor(ind1)+indf;
                
                if debug
                    figure(1); hold on
                    plot(xi,yi,'r*')
                end
                
                
                % If hit a corner, we will consider both edges
                if ~isempty(corner)
                    edge = corner;
                end
                
                % Store this edge
                if ~lastEdgeReset
                    lastedge = edge;
                else
                    lastEdgeReset = 0;
                end
                lastcolr = lastcol;
                lastcol = 0;
                justreflect = 1;
  
                
                %if an edge is hit, check what class of frame it comes from
                switch class(frmout)
                    case 'causticsFrame'
                        
                        switch class(frmout)
                            case 'causticsFrame'
                                
                                if ~all(frmout.edgestyle(edge))
                                    % Case of non-ohmic edge
                                    
                                    if rand<p_ifbounce_then_scatter
                                        % Scatter
                                        ind2 = BS.inject(frmout.edgenorm(edge));
                                        
                                        % Store coordinates
                                        % Before Scattering
                                        ind1s(2,k-1) = indi;
                                        % After Scattering
                                        x1s(k) = xi;
                                        y1s(k) = yi;
                                        ind1s(1,k) = ind2;
                                        
                                        % Update coordinates
                                        x2 = xi; y2 = yi;
                                        
                                    else
                                        % Specular reflection
                                        edgenorm = frmout.edgenorm(edge);
                                        %inprobEdge = find(BS.unorms == edgenorm);
                                        
                                        [~,~,ind2,~] = BS.specReflect(x1,y1,xi,yi,indi,frmout.edgenorm(edge));
                                        if isnan(ind2)
                                            fprintf('ind2 is Nan\n')
                                        end
                                        
                                        % Store coordinates
                                        % Before reflection
                                        ind1s(2,k-1) = indi;
                                        % After reflection
                                        x1s(k) = xi;
                                        y1s(k) = yi;
                                        ind1s(1,k) = ind2;
                                        
                                        % Update coordinates
                                        x2 = xi; y2 = yi;
%                                         
%                                         if BS.fracInd(BS.inprobCheck(inprobEdge,:),indo) ~= 0
%                                             % Store coordinates
%                                             % Before reflection
%                                             ind1s(2,k-1) = indi;
%                                             % After reflection
%                                             x1s(k) = xi;
%                                             y1s(k) = yi;
%                                             ind1s(1,k) = indo;
%                                             
%                                             % Update coordinates
%                                             x2 = xi; y2 = yi; ind2 = indo;
%                                         else
%                                             fprintf('Bad reflection occured\n')
%                                             fprintf('Steps since last collison: %d\n',lastcolr)
%                                             k = k-1;
%                                             if frm1.inbounds(x2,y2) % Check if next step in
%                                                 fprintf('Next step was inbounds, continuing\n')
%                                                 lastedge = 0;
%                                             else
%                                                 fprintf('Fatal error, reinjecting\n')
%                                                 error = 1; % Set error flag
%                                                 active_trajectory=0;
%                                                 
%                                                 k = kstore; % Reset k to value before current electron
%                                                 bad_x1s = x1s;
%                                                 bad_y1s = y1s;
%                                                 bad_ind1s = ind1s;
%                                                 bad_BS = BS;
%                                                 
%                                                 bad_arcdata = [bad_x1s' bad_y1s' bad_ind1s'];
%                                                 lastind = find(bad_arcdata(:,3) ~= 0,1,'last');
%                                                 bad_arcdata(lastind+1:end,:)=[];
%                                             end
%                                     end
                                    
                                        
                                        x1r = x1; y1r = y1; ind1r = ind1;
                                        xir = xi; yir = yi; indir = indi;
                                        ind2r = ind2;
                                    end
                                k=k+1;
                                
                                elseif frmout.terminal_ohmic_test(edge)
                                    %elseif sum(frmout.edgestyle(edge)==terminal_ohmics)
                                    % Case of grounded ohmic edge
                                    
                                    if rand<p_bounce_off_ohmic
                                        % Bounces off ohmic
                                        
                                        if rand<p_ifbounce_then_scatter
                                            % Scatter
                                            ind2 = BS.inject(frmout.edgenorm(edge));
                                            
                                            % Store coordinates
                                            % Before Scattering
                                            ind1s(2,k-1) = indi;
                                            % After Scattering
                                            x1s(k) = xi;
                                            y1s(k) = yi;
                                            ind1s(1,k) = ind2;
                                            
                                            % Update coordinates
                                            x2 = xi; y2 = yi;
                                            
                                        else
                                            % Specular reflection
                                            [~,~,ind2,~] = BS.specReflect(x1,y1,xi,yi,indi,frmout.edgenorm(edge));
                                            

                                            
                                            % Store coordinates
                                            % Before reflection
                                            ind1s(2,k-1) = indi;
                                            % After reflection
                                            x1s(k) = xi;
                                            y1s(k) = yi;
                                            ind1s(1,k) = ind2;
                                            
                                            % Update coordinates
                                            x2 = xi; y2 = yi;
                                        end
                                        k=k+1;
                                        
                                    else
                                        % Sunk to ground
                                        active_trajectory=0;
                                        
                                        % Pick one of the terminal ohmic
                                        % edges if corner
                                        edge = edge(find(frm1.edgestyle(edge)~=0,1));
                                        
                                        % Store final index
                                        ind1s(2,k-1) = indi;
                                    end
                                    
                                else
                                    % Case of floating ohmic edge
                                    
                                    if rand<p_bounce_off_ohmic
                                        % Bounce of ohmic
                                        
                                        if rand<p_ifbounce_then_scatter
                                            % Scatter
                                            ind2 = BS.inject(frmout.edgenorm(edge));
                                            
                                            % Store coordinates
                                            % Before Scattering
                                            ind1s(2,k-1) = indi;
                                            % After Scattering
                                            x1s(k) = xi;
                                            y1s(k) = yi;
                                            ind1s(1,k) = ind2;
                                            
                                            % Update coordinates
                                            x2 = xi; y2 = yi;
                                            
                                        else
                                            % Specular reflection
                                            [~,~,ind2,~] = BS.specReflect(x1,y1,xi,yi,indi,frmout.edgenorm(edge));
                                            

                                            
                                            % Store coordinates
                                            % Before reflection
                                            ind1s(2,k-1) = indi;
                                            % After reflection
                                            x1s(k) = xi;
                                            y1s(k) = yi;
                                            ind1s(1,k) = ind2;
                                            
                                            % Update coordinates
                                            x2 = xi; y2 = yi;
                                        end
                                        k=k+1;
                                        
                                    else
                                        % Absorbs and reinjects
                                        repedge = edge(find(frm1.edgestyle(edge)~=0,1));
                                                                              
                                        edgenum(repedge)=edgenum(repedge)+1;
                                        [x2,y2,norm] = frmout.ogs{frm1.edgestyle(repedge)}.inject_position();
                                        edge = frm1.ogs{frm1.edgestyle(repedge)}.nedge_in_ohmic(frm1.ogs{frm1.edgestyle(repedge)}.edge_inject);
                                        lastedge=edge; % Sets last edge to the injector edge
                                        ind2 = BS.inject(norm);
                                        % Store coordinates
                                        % Before absorption
                                        ind1s(2,k-1) = indi;
                                        % Reinjection
                                        x1s(k) = x2;
                                        y1s(k) = y2;
                                        ind1s(1,k) = ind2;
                                        k=k+1;
                                    end
                                end
                                
                        end
                        %                 case 'causticsBarrier'
                        %                     %refracts or reflects off of barrier
                        %                     [theta1,B_cur]=frmout.reflect_or_transmit(arc1(k).incident_angle(x1,y1),B_cur0,edge);
                        %                     k=k+1;
                        
                end
            end % end collision loop
        catch
            fprintf('An error occured during the current electron\n')
            %display(sstore)
            %display(store)
            error = 1; % Set error flag
            active_trajectory=0;
            
            k = kstore; % Reset k to value before current electron
            bad_x1s = x1s;
            bad_y1s = y1s;
            bad_ind1s = ind1s;
            bad_BS = BS;
            
            bad_arcdata = [bad_x1s' bad_y1s' bad_ind1s'];
            lastind = find(bad_arcdata(:,3) ~= 0,1,'last');
            bad_arcdata(lastind+1:end,:)=[];
            
        end
    end % While active_trajectory
    
    %k=k+1;
    
    if ~error
        edgenum(edge)=edgenum(edge)+1;
        
        if(~sum(edge==frm1.ogs{injector_ohmic}.nedge_in_ohmic))
            n_in=n_in+1;
        end
    end
end


% Generate final form of arcdata that will be returned
% Narcs=length(arc1);

arcdata = [x1s' y1s' ind1s'];
lastind = find(arcdata(:,3) ~= 0,1,'last');
arcdata(lastind+1:end,:)=[];



