Xmin=-3.5;
Xmax=3.5;
Ymin=-90;
Ymax=90;
DeltaX=.01;
DeltaY=DeltaX;
% Step size along arcs (should be smaller than DeltaX)
dL=.01;
% Color contrast for viewing
Contrast=10;

% Set to true to export png files
save_plots=true;
% Frame rate of gif output
frame_rate=5;

%%
% Computes the meshgrid data for the colormap plot




addpath([savepath,fname])
load([savepath,fname,'/frameData.mat']);

%find the files associated with simulation results

filelist = ls([savepath,fname,'/*arcdata.txt']);

classlist = ls([savepath,fname,'/*classes.mat']);

Nfiles = size(filelist,1);% number of files to loop through


myCmaps;

interp = 0;
scaling = 0.1;

%generate X,Y meshgrid for binning arc points
[X,Y]=meshgrid(Xmin:DeltaX:Xmax,Ymin:DeltaY:Ymax);
save([savepath,fname,'/imageProperties.mat'],'X','Y','dL');

%set up arrays and offsets for rounding coordinates to nearest grid point
X1D=squeeze(X(1,:));
Y1D=squeeze(Y(:,1));
Xoff=find(abs(X1D)==min(abs(X1D)));
Yoff=find(abs(Y1D)==min(abs(Y1D)));

 

%for k = 585
for k = 1:Nfiles
    
    t=tic; %timing for loop iteration
    
    % Load in data file
    fid=fopen([savepath,fname,'/',squeeze(filelist(k,:))]); % open sim. file
    arcdata=fread(fid,'double');
    arcdata=reshape(arcdata,[],4);
    fclose(fid);
    
    % Load in arcclass
    load([savepath,fname,'/',squeeze(classlist(k,:))]);
    
    %initial Z matrix
    Z=zeros(size(X));
    
    Narcs=size(arcdata,1);
    
    %Step through list of arcs.
    for i=1:Narcs
        
        [xout,yout] = arcclass.arc(arcdata(i,:));
        
        if interp
            % Initialize interpolated arrays
            xouti = xout(1);
            youti = yout(1);
            
            Del = scaling*min([DeltaX, DeltaY]);
            
            for j = 2:length(xout)
                Nx = ceil(abs(xout(j)-xout(j-1))./Del);
                Ny = ceil(abs(yout(j)-yout(j-1))./Del);
                N = [Nx;Ny];
                [Nmax, ~] = max(N(:,end));
                
                newX = linspace(xout(j-1),xout(j),Nmax);
                newY = linspace(yout(j-1),yout(j),Nmax);
                xouti = [xouti newX(2:end)];
                youti = [youti newY(2:end)];
            end
            
            xout = xouti;
            yout = youti;
        end
        
        %step throug the list of coordinates and round to nearest grid
        %point
        for n=1:length(xout)
            if(xout(n)>=Xmin&&xout(n)<=Xmax&&yout(n)>=Ymin&&yout(n)<=Ymax)
                npoint=round(yout(n)/DeltaY+Yoff);
                mpoint=round(xout(n)/DeltaX+Xoff);
                
                %increment a value of 1 for the binned location of the arc
                Z(npoint,mpoint)=Z(npoint,mpoint)+1;
            end
        end
        
    end
    %save Z matrix for later viewing
    save([savepath,fname,'/',fname,'_',num2str(k,'%05d'),'Z.mat'],'Z');
    
    %display image
    image(X(1,:),Y(:,1),min(Z,N_inject/Contrast)/(N_inject/Contrast)*255);
    set(gca,'ydir','normal');
    axis equal;
    hold on;
    shading flat;
    
    colormap(cmap.red); %optional colormap from "myCmaps.m"
    %set(gca,'clim',[0 N_inject/Contrast]);
    
    %plot the frame and ohmics for visual reference
    plotFrameGroup(frmgrp);
    
    hold off
    drawnow();
    
    t2=toc(t);
    fprintf('%d/%d:  %0.1f\n',k,Nfiles,t2); %output the computation stats
end

