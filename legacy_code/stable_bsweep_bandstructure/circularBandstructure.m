classdef circularBandstructure < handle
    
    % b ~ 0.052 Angstrom^-1
    % vb = 340 meV
    
    
    properties
        
        % Fitting parameters of Fermi surface
        k0  = 1;  % Å-1
        vf  = 4.8;   % eVÅ (7.3*10^5 m/s),
        
        % Input parameters
        phi = 0;     % Moire angle wrt sample
        
        % Computed values
        K            % K vectors of trajectory
        K_center     % Center of K space trajectory
        K_n          %
        v_g          % Group velocity at each K point
        R            % Radial position of real space trajecotry
        theta        % theta of realspace trajectory
        Rx           % X position of trajectory
        Ry           % Y position of trajectory
        dRx          % Step to next Rx
        dRy          % Step to next Ry
        unorms       % Unique norms in geometry
        inprob       % Injection prob for all unique norms
        inprobCheck  % Used for checking prob of any index
        B            % Used for switching direction of trajectory
        
        % Fermi Surface generation
        n_theta = 1000; % Number of points for generating Fermi surface
        n_interp = 1000; % Number of points to interpolate Fermi surface
        
    end %properties
    
    
    methods
        %% Load in and reshape mesh
        function init(obj)
            
            
            % Calculate Fermi surface
            if obj.B>0
                theta_loop = linspace(0,2*pi,obj.n_theta);
            else
                theta_loop = linspace(2*pi,0,obj.n_theta);
            end
            kf=obj.k0*ones(size(theta_loop));
            
            K_loop = [kf.*cos(theta_loop); kf.*sin(theta_loop)]; % To Cartesian
            
            approxMid = mean(K_loop,2);
            
            K_loop = obj.rotate(K_loop,obj.phi); % Rotate xtal axis
                        
            
            arcL = sqrt(sum(diff(K_loop,[],2).^2,1)); % Arc length between points % Arc length between points
            arcL = [0 arcL]; % Append starting point
            cumArcL = cumsum(arcL);
            interpPoints = linspace(0,cumArcL(end),obj.n_interp); % Points for interpolation
            obj.K = interp1(cumArcL,K_loop',interpPoints,'linear')';
            
            
            obj.K_n = size(obj.K,2);
            obj.K_center = mean(obj.K,2);
            finalArcL =  sqrt(sum(diff(obj.K,[],2).^2,1));
   
        end
        

        
        %% Calculate the orbit in real space by scaling and rotating the orbit
        function KtoR(obj,B)
            % Takes in magnitude of field
            % It is already expecting a K orbit computed using
            % the correct sign field
            
            % Field in units of T*C/hbar
            
            tmpr = (1/B).*(obj.K-repmat(obj.K(:,1),1,size(obj.K,2)));
            obj.Rx = tmpr(2,:) - mean(tmpr(2,:));
            obj.Ry = -tmpr(1,:) + mean(tmpr(1,:));
            %             if B<0
            %                 obj.Rx = fliplr(obj.Rx);
            %                 obj.Ry = fliplr(obj.Ry);
            %                 obj.v_g = -fliplr(obj.v_gstore);
            %             else
            %                 obj.v_g = obj.v_gstore;
            %             end
            
            obj.dRx = obj.Rx(2:end) - obj.Rx(1:end-1);
            obj.dRy = obj.Ry(2:end) - obj.Ry(1:end-1);
            obj.R = sqrt(obj.Rx.^2+obj.Ry.^2);
            obj.theta = atan2(obj.Ry,obj.Rx);
            
            obj.Rx = round(obj.Rx,4);
            obj.Ry = round(obj.Ry,4);
            obj.dRx = round(obj.dRx,4);
            obj.dRy = round(obj.dRy,4);
            
        end
        

        
        %% Plot K and Real space trajectories
        function plotOrbits(obj, fignum)
            % Plot K and R Space Trajectories
            figure(fignum); clf;
            subplot(1, 2, 1)
            plot(obj.K(1,:),obj.K(2,:),'.-')
            xlabel('$K_x(t) (\AA^{-1})$','interpreter','latex','fontsize',18)
            ylabel('$K_y(t) (\AA^{-1})$','interpreter','latex','fontsize',18)
            grid on
            axis tight
            
            subplot(1, 2, 2)
            hold on
            plot(obj.Rx,obj.Ry,'.-')
            xlabel('$x(t) (\mu m)$','interpreter','latex','fontsize',18)
            ylabel('$y(t) (\mu m)$','interpreter','latex','fontsize',18)
            grid on
            axis tight
            
        end
        
        
        %% Returns arc for given data
        function [xout,yout] = arc(obj,arcdata)
            x1 = arcdata(1);
            y1 = arcdata(2);
            ind1 = arcdata(3);
            ind2 = arcdata(4);
            

            
            xout = x1;
            yout = y1;
            i = ind1; % Initialize index
            n = 2; % Indexing for output
            
            if 1 > (ind1-ind2) &&  (ind1-ind2) >= 0
                xout = x1;
                yout = y1;
            else
                % If starting index is fractional
                if ceil(i) ~= i && floor(i) ~= floor(ind2)
                    xout(n) = xout(n-1) + obj.Rx(ceil(i)) - obj.fracInd(obj.Rx,i);
                    yout(n) = yout(n-1) + obj.Ry(ceil(i)) - obj.fracInd(obj.Ry,i);
                    indtest(n-1) = i;
                    i = mod(floor(i),size(obj.Rx,2)-1)+1;
                    n = n + 1;
                end
                
                while i ~= floor(ind2) && i == round(i)
                    inext = mod(i,size(obj.Rx,2)-1)+1;
                    xout(n) = xout(n-1)+obj.Rx(inext)-obj.Rx(i);
                    yout(n) = yout(n-1)+obj.Ry(inext)-obj.Ry(i);
                    %indtest(n-1) = i;
                    i = inext;
                    n = n + 1;
                    if ind2 == size(obj.Rx,2)
                        ind2 =1;
                    end
                end
                %i = mod(i-2,size(obj.dRx,2))+1;
                
                
                if i ~= ind2
                    xout(n) = xout(n-1)+obj.fracInd(obj.Rx,ind2)-obj.fracInd(obj.Rx,i);
                    yout(n) = yout(n-1)+obj.fracInd(obj.Ry,ind2)-obj.fracInd(obj.Ry,i);
                    %indtest(n-1) = ind2;
                    %                i = ceil(i);
                    %                n = n + 1;
                end
            end
            
            
            
            
        end
        
        %% Compute the probability distribution for each unique edge
        function injectProb(obj,edgenorms)
            obj.unorms = unique(edgenorms);
            
            S = sqrt(obj.dRx.^2 + obj.dRy.^2);
            Smax = max(S);
            ang = atan2(obj.dRy,obj.dRx);
            
            for i = 1:length(obj.unorms)
                % probtmp = cos(ang-obj.unorms(i)) .* S/Smax;
                probtmp = cos(ang-obj.unorms(i));
                probtmp(probtmp<0) = 0;
                obj.inprob(i,:) = probtmp;
            end
            
            obj.inprobCheck = [obj.inprob obj.inprob(:,1)];
            
        end
        
        %% Injection code
        function n = inject(obj,edgenorm)
            %             ind = find(edgenorm == obj.unorms);
            %             r = rand;
            %             n = find(r <= obj.inprob(ind,:),1);
            n = [];
            
            if length(edgenorm) == 1
                edge = find(edgenorm == obj.unorms);
                while isempty(n)
                    ind = randi([1,size(obj.inprob,2)-1]);
                    r = rand;
                    if r < obj.inprob(edge,ind)
                        n = ind;
                    end
                end
            else
                edge1 = find(edgenorm(1) == obj.unorms);
                edge2 = find(edgenorm(2) == obj.unorms);
                while isempty(n)
                    ind = randi([1,size(obj.inprob,2)-1]);
                    r = rand;
                    if r < obj.inprob(edge1,ind)
                        n = ind;
                    end
                    if obj.inprob(edge1,n) == 0
                        n = [];
                    end
                    
                end
            end
            
        end
        

        
        
        %% Rotates a given input vector by ang
        function out = rotate(~,in,ang)
            % Rotates coordinate pairs by theta (positive is CCW)
            out(1,:)=cos(ang)*in(1,:)-sin(ang)*in(2,:);
            out(2,:)=cos(ang)*in(2,:)+sin(ang)*in(1,:);
            
        end
        
        
         %% Integrates phase propertly to remove mod(2*pi)
        function out = phaseunwind(~,in)
            
            N=length(in);
            
            del=in(2:end)-in(1:(end-1));
            
            del(find(del>2*pi*.8))=100; %#ok<*FNDSB>
            del(find(del<-2*pi*.8))=-100;
            del(find(abs(del)<90))=0;
            del=[0;del];
            
            out=in-cumsum(del*2*pi/100);
            
        end
        
        %% finds interpolated zerocrossing
        function out = zeroCrossInterp(~,varargin)
            
            if nargin==3
                X=varargin{1};
                Y=varargin{2};
                
            end
            
            if nargin==2
                Y=varargin{1};
                X=(1:length(Y))';
            end
            
            zerocrosses=find(Y(2:end).*Y(1:end-1)<0);
            zerocrosses = [zerocrosses; find(Y == 0)];
            
            if isempty(zerocrosses)
                out = [];
            else
                for i = 1:length(zerocrosses)
                    ind = zerocrosses(i);
                    if ind == X(end)
                        % Only happens when 0 exactly at end
                        xout(i) = X(ind);
                    elseif Y(ind) == Y(ind+1) 
                        % If two adjacent exact zeros, don't want NAN
                        xout(i) = X(ind);
                    else
                        xout(i)=X(ind)-Y(ind).*(X(ind)-X(ind+1))./(Y(ind)-Y(ind+1));
                    end
                end
                %exactcrosses=X(Y==0);
                %out=sort([xout;exactcrosses]);
                out = sort(xout);
            end
            
        end
        
        
        %% Finds intersection with a boundary
        function [xi,yi,indi] = intersect(~,x2,x1,y2,y1,px2,px1,py2,py1,ind1)
            % This all assumes linear interpolation between points
            % We are always stepping from ind1 to an integer and need
            % to find indi
            [xi,yi] = polyxpoly([x2 x1],[y2 y1],[px2 px1],[py2,py1]);
            ind1 = mod(ind1,1);
            if x2 ~= x1
                indi = ind1 + (1-ind1)*(xi-x1)/(x2-x1);
            else
                indi = ind1 + (1-ind1)*(yi-y1)/(y2-y1);
            end
            
            
            % In case it intersects the same edge twice
            % Take the point closest to the original
            [indi, i] = min(indi);
            xi = xi(i);
            yi = yi(i);
            
        end
        
                %% Finds intersection with a boundary
        function [indi] = intersectInd(~,x1,xi,x2,y1,yi,y2,ind1)
            % This all assumes linear interpolation between points
            % We are always stepping from ind1 to an integer and need
            % to find indi
            ind1 = mod(ind1,1);
            if x2 ~= x1
                indi = ind1 + (1-ind1)*(xi-x1)/(x2-x1);
            else
                indi = ind1 + (1-ind1)*(yi-y1)/(y2-y1);
            end
            
            
            % In case it intersects the same edge twice
            % Take the point closest to the original
            [indi, i] = min(indi);
            xi = xi(i);
            yi = yi(i);
            
        end
        
        
        %% Specular reflection
        function [x2,y2,no,n2] = specReflect(obj,x1,y1,xi,yi,ni,edgenorm)
            
            
            edgenormstore = edgenorm;
            nedge = length(edgenorm);
            meanNormFlag = 0;
            if nedge>1
                ang = atan2(yi-y1,xi-x1);
                test = abs(ang-edgenorm)>pi/2;
                if test(1)
                    edgenorm = edgenorm(1);
                elseif test(2)
                    edgenorm = edgenorm(2);
                elseif all(test)
                    % Bounce off both edges
                    edgenorm = mean(edgenorm);
                    meanNormFlag = 1;
                elseif edgenorm(1) == edgenorm(2)
                    % hitting a boundary of edgestyles
                    edgenorm = edgenorm(1);
                else
                    fprintf('Fatal error in specReflect edgenorm\n')
                end
            end
            
            [Rxo,Ryo,nr] = obj.reflectstates(ni,edgenorm);
            
            nround = round(nr);
            
            % Want to selected the closest state with finite injection
            % probability
            
            try
                if ~meanNormFlag % only have checks for reflecting off of single edge
                    edge = find(edgenorm == obj.unorms);
                    i = 1;
                    no = [];
                    while isempty(no)
                        if obj.inprobCheck(edge,nround(i))~=0
                            no = nround(i);
                            Rxo = Rxo(i);
                            Ryo = Ryo(i);
                        end
                        i = i+1;
                    end
                else % Bouncing off of two edges
                    edge1 = find(edgenormstore(1) == obj.unorms);
                    edge2 = find(edgenormstore(2) == obj.unorms);
                    i = 1;
                    no = [];
                    while isempty(no)
                        if obj.inprobCheck(edge1,nround(i))~=0 && obj.inprobCheck(edge2,nround(i))~=0
                            no = nround(i);
                            Rxo = Rxo(i);
                            Ryo = Ryo(i);
                        end
                        i = i+1;
                    end
                end
            catch
                % We will expand on the above search
                nf = floor(nr);
                nc = ceil(nr);
                ncheck = zeros(length(nf) + length(nc),1);
                ncheck(1:2:end) = nf; ncheck(2:2:end) = nc;
                if ~meanNormFlag % only have checks for reflecting off of single edge
                    edge = find(edgenorm == obj.unorms);
                    i = 1;
                    no = [];
                    while isempty(no)
                        if obj.inprobCheck(edge,ncheck(i))~=0
                            no = ncheck(i);
                            Rxo = [];
                            Ryo = [];
                        end
                        i = i+1;
                    end
                else % Bouncing off of two edges
                    edge1 = find(edgenormstore(1) == obj.unorms);
                    edge2 = find(edgenormstore(2) == obj.unorms);
                    i = 1;
                    no = [];
                    while isempty(no)
                        if obj.inprobCheck(edge1,ncheck(i))~=0 && obj.inprobCheck(edge2,ncheck(i))~=0
                            no = ncheck(i);
                            Rxo = [];
                            Ryo = [];
                        end
                        i = i+1;
                    end
                end

            end
            
            n2 = ceil(no);
            x2 = xi + obj.Rx(n2) - Rxo;
            y2 = yi + obj.Ry(n2) - Ryo;


        end
        
        %% Reflects a given state into the closest state in K space
        function [Rxout,Ryout,nout] = reflectstates(obj,ni,edgenorm)
            % Returns all possible reflected states sorted by distance from
            % orignal point
            
            
            % We will reflect parallel to the wall in our real space orbit
            %m = (y1-y0)/(x1-x0); % Slope of the wall
            %gamma = atan(m);
            gamma = edgenorm + pi/2;
            
            R_rot = obj.rotate([obj.Rx; obj.Ry],-gamma);
            b = (1-mod(ni,1))*R_rot(2,floor(ni)) + mod(ni,1)*R_rot(2,ceil(ni));
            
            ncross = obj.zeroCrossInterp((R_rot(2,:) - b)'); % Interpolate the zero crossings
            ncross = sort(ncross);
            if ncross(1) == 1 && ncross(end) == length(obj.Rx)
               ncross(end) = []; 
            end
            
            
            % Weighted average to find exact position of reflected state

            for i = 1:length(ncross)
                Rxcross(i) = obj.fracInd(obj.Rx,ncross(i));
                Rycross(i) = obj.fracInd(obj.Ry,ncross(i));
                
            end

            
            Rx_n = obj.fracInd(obj.Rx,ni);
            Ry_n = obj.fracInd(obj.Ry,ni);
            
            Rcross = sqrt(Rxcross.^2+Rycross.^2);
            thetacross = atan2(Rycross,Rxcross);
            d = sqrt((Rx_n - Rxcross).^2 + (Ry_n - Rycross).^2);
            %             dthresh = 1E-8;
            %             if isempty(min(d(d>eps)))
            %                 fprintf('Error: Did not find a reflected state within threshold\n')
            %                 fprintf('%f\n',d)
            %             end
            [~, ind] = sort(d);
            
           
            
            
            Rxout = Rxcross(ind);
            Ryout = Rycross(ind);
            Rout = Rcross(ind);
            thetaout = thetacross(ind);
            nout = ncross(ind);
            
        end
        
        %% Linear interpolation to find value at a fractional index
        function out = fracInd(~,in,n)
            out = (1-mod(n,1))*in(floor(n)) + ...
                mod(n,1)*in(ceil(n));
        end
        
        %% Returns a colorplot of the Injection Probabilities
        function plotInProb(obj,edgenum,fignum)
            figure(fignum); clf;
            plot(obj.inprob(edgenum,:))
            
            figure(fignum+1); clf;
            subplot(1, 2, 1)
            sz = size(obj.inprob(edgenum,:));
            scatter3(obj.K(1,1:end-1),obj.K(2,1:end-1),zeros(sz),15*ones(sz),obj.inprob(edgenum,:),'filled')
            xlabel('$K_x$','interpreter','latex','fontsize',18)
            ylabel('$K_y$','interpreter','latex','fontsize',18)
            colormap(jet)
            colorbar
            view(0,90)
            axis tight
            grid on
            
            subplot(1, 2, 2)
            hold on
            scatter3(obj.Rx(1:end-1),obj.Ry(1:end-1),zeros(sz),15*ones(sz),obj.inprob(edgenum,:),'filled')
            xlabel('$x(t)$','interpreter','latex','fontsize',18)
            ylabel('$y(t)$','interpreter','latex','fontsize',18)
            colorbar
            view(0,90)
            axis tight
            grid on
            
        end
        
    end %methods
    
    
end








































