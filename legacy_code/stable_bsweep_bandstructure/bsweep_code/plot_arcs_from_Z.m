%load properties of the simulation and meshgrids
load([savepath,fname,'/frameData.mat']);
load([savepath,fname,'/imageProperties.mat']);
minX=min(min(X)); minY=min(min(Y));
maxX=max(max(X)); maxY=max(max(Y));

%find the available processed Z matricies
addpath([savepath,fname])
filelist=ls([savepath,fname,'/*Z.mat']);
Nfiles=size(filelist,1);%number of files to loop through

%load in my homemade cmaps choices are teal and red
myCmaps;


figure(1);
h=subplot(1,1,1);
if save_plots
    mkdir([savepath, fname, '/PNG'])
end


for i=1:Nfiles
    
    cla(h);
    
    hold off;
    %load processed Z matrix
    load([savepath,fname,'/',fname,'_',num2str(i,'%05d'),'Z.mat']);
    
    %using image. Z is thresholded and normalized to 255
    image(X(1,:),Y(:,1),min(Z,N_inject/Contrast)/(N_inject/Contrast)*255);
    set(gca,'ydir','normal')
    hold on
    
    %plot the frame and ohmics for visual reference
    plotFrameGroup(frmgrp);
    
    %use cmap of your choice
    colormap(cmap.red);
    
    drawnow();
    if save_plots
        
        im = frame2im(getframe(gcf()));
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind,cm,[savepath, fname, '/PNG/',fname,'_', num2str(i), '_PNG.png'])  
        if i == 1
            imwrite(imind,cm,[savepath, fname, '/',fname,'.gif'],'gif', 'Loopcount',Inf,'DelayTime',1/frame_rate);
        else
            imwrite(imind,cm,[savepath, fname, '/',fname,'.gif'],'gif','WriteMode','append','DelayTime',1/frame_rate);
        end
        
    end
end

