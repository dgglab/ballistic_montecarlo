addpath('..')
% dirpath = '..\..\Data\stable_bsweep_bandstructure\miniBand\';
% savedir = 'phi_0 B_min_-0.0075959 B_max_0.0075959';
% savepath = [dirpath,savedir,'\xip1\'];
% dirlist = ls([dirpath,savedir,'\xip1\phi*']);
dirpath = '..\..\Data\stable_bsweep_bandstructure\delafossite\Sherlock\';
% savedir = 'Ef=m188_scattering';
% savepath = [dirpath,savedir,'\xim1Data\'];
% dirlist = ls([dirpath,savedir,'\xim1Data\phi*']);
%savedir = 'phi_0 B_min_-3 B_max_0 B_N_1001 N_301 p_scatter_1 p_ohmic_0_run2\';
savedir = 'rotData\';
savepath = [dirpath,savedir];
dirlist = ls([dirpath,savedir,'\phi*']);
Ndir = size(dirlist,1);

dir = strtrim(squeeze(dirlist(1,:)));

fname = dir;

%load properties of the simulation and meshgrids
load([savepath,fname,'/frameData.mat']);
load([savepath,fname,'/imagePropertiesTest.mat']);
minX=min(min(X)); minY=min(min(Y));
maxX=max(max(X)); maxY=max(max(Y));

% Set to true to export png files
save_plots=true;
% Frame rate of gif output
frame_rate=5;

Contrast=1.5;
skip = 1;

% for iDir = 1:Ndir
%Ndir=1;
for iDir = 1:Ndir
    dir = strtrim(squeeze(dirlist(iDir,:)));
    display(dir)
    
    dirname = dir;
    
    
    %find the available processed Z matricies
    addpath([savepath,dirname])
    filelist=ls([savepath,dirname,'/*Zpub.mat']);
    
    Nfiles=size(filelist,1);% number of files to loop through
    
    %load in my homemade cmaps choices are teal and red
    myCmaps;
    
    
    figure(1);
    h=subplot(1,1,1);
    
    
    fname = strtrim(squeeze(filelist(1,:)));
    
    load([savepath,dirname,'/',fname]);
   
    
    Ztemp = zeros(size(Z,1),size(Z,2),length(1:skip:Nfiles));
    if iDir == 1
        Zstore = Ztemp;
    end
    
    j = 1;
    
    %for i=1:skip:Nfiles
    for i=1
        fname = strtrim(squeeze(filelist(i,:)));
        
        cla(h);
        
        hold off;
        %load processed Z matrix
        load([savepath,dirname,'/',fname]);
        Ztemp(:,:,j) = Z;
        j = j + 1;
        
    end
    
    
    Zstore = Zstore + Ztemp/Ndir;
end

%%
if save_plots
    mkdir([savepath, '/PNG'])
    mkdir([savepath, '/GIF'])
end
dirname = strtrim(squeeze(dirlist(1,:)));
for i = 1:size(Zstore,3)
    Zstore = imgaussfilt(Zstore(:,:,i),3);
    %using image. Z is thresholded and normalized to 255
    image(X(1,:),Y(:,1),(min(Zstore(:,:,i),N_inject/Contrast)/(N_inject/Contrast)*255));
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
        imwrite(imind,cm,[savepath, '/PNG/',dirname,'_', num2str(i), '_PNG.png'])
        if i == 1
            imwrite(imind,cm,[savepath, '/GIF/',dirname,'.gif'],'gif', 'Loopcount',Inf,'DelayTime',1/frame_rate);
        else
            imwrite(imind,cm,[savepath, '/GIF/',dirname,'.gif'],'gif','WriteMode','append','DelayTime',1/frame_rate);
        end
        
    end
end