% Author - Aaron Sharpe (email: asharpe1@stanford.edu)
% Description - This function creates a color map of two different data
% sets acording to a 2D colormap


% Inputs
skip = 1; % Manually skip frames

fixPlot = 1; % Manually removes points outside the boundary

g = 1/3; %gamma power factor for contrast in colormap

save_plots=true; % Set to true to export png files

frame_rate=5; % Frame rate of gif output


addpath('..')
dirpath = '..\..\Data\stable_bsweep_bandstructure\miniBand\';
savedir = 'Ef=m188';
subdir = '\xip1Data\';
savepath = [dirpath,savedir,subdir];
dirlist = ls([dirpath,savedir,subdir,'phi*']);
Ndir = size(dirlist,1);

dir = strtrim(squeeze(dirlist(1,:)));

addpath('..')
dirpath2 = '..\..\Data\stable_bsweep_bandstructure\miniBand\';
savedir2 = 'Ef=m188';
subdir2 = '\xim1Data\';
savepath2 = [dirpath2,savedir2,subdir2];
dirlist2 = ls([dirpath2,savedir2,subdir2,'phi*']);
Ndir2 = size(dirlist2,1);

%% Load in some image properties

dir = strtrim(squeeze(dirlist(1,:)));

% Load properties of the simulation and meshgrids
load([savepath,dir,'/frameData.mat']);
load([savepath,dir,'/imagePropertiesTest.mat']);
minX=min(min(X)); minY=min(min(Y));
maxX=max(max(X)); maxY=max(max(Y));

%% Load in first directory
for iDir = 1:Ndir
    dir = strtrim(squeeze(dirlist(iDir,:)));
    display(dir)
    
    dirname = dir;
    
    
    %find the available processed Z matricies
    addpath([savepath,dirname])
    filelist=ls([savepath,dirname,'/*Ztest.mat']);
    
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
    
    for i=1:skip:Nfiles
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

Zstore1 = Zstore;

%% Load in second directory

dirpath = dirpath2;
savedir = savedir2;
subdir = subdir2;
savepath = savepath2;
dirlist = dirlist2;
Ndir = Ndir2;

dir = strtrim(squeeze(dirlist(1,:)));

%load properties of the simulation and meshgrids
load([savepath,dir,'/frameData.mat']);
load([savepath,dir,'/imagePropertiesTest.mat']);
minX=min(min(X)); minY=min(min(Y));
maxX=max(max(X)); maxY=max(max(Y));



for iDir = 1:Ndir
    dir = strtrim(squeeze(dirlist(iDir,:)));
    display(dir)
    
    dirname = dir;
    
    
    %find the available processed Z matricies
    addpath([savepath,dirname])
    filelist=ls([savepath,dirname,'/*Ztest.mat']);
    
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
    
    for i=1:skip:Nfiles
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

Zstore2 = Zstore;


%% Manually set all pixels outside the device body to zero for the electrons that escaped :(
if fixPlot
    in = inpolygon(X,Y,frmgrp.px,frmgrp.py);
    for i = 1:size(Zstore1,3)
        Zstore1(:,:,i) = Zstore1(:,:,i).*in;
        Zstore2(:,:,i) = Zstore2(:,:,i).*in;
    end
end


%%
if save_plots
    mkdir([dirpath,savedir, '/PNG'])
    mkdir([dirpath,savedir, '/GIF'])
end

% Create the Color map for later
AminBmin = [1, 1, 1]; % No count in either channel
AmaxBmin = [0.8, 0, 0]; % Max count in channel A
AminBmax = [0, 0, 0.8]; % Max count in channel B
AmaxBmax = [0.4, 0, 0.4]; % Max count in both channels
CMAP  = Colormap2d(500, 500, AminBmin, AminBmax, AmaxBmin, AmaxBmax);
%figure; imshow(CMAP)


% Adjust Range of input data to be inbetween 0 and 1
Zmin = min([min(Zstore1(:)),min(Zstore2(:))]);
Z1 = Zstore1 - Zmin;
Z2 = Zstore2 - Zmin;
Zmax = max([max(Zstore1(:)),max(Zstore2(:))]);
Z1 = Z1/Zmax;
Z2 = Z2/Zmax;

%Adjust contrast for better viewing, for quick good looking images using 
%non consistent color scale, uncomment imadjust later

Z1 = Z1.^g;
Z2 = Z2.^g;


for i = 1:size(Zstore,3)
    
    %Remove intensity map for time i
    Z1s = Z1(:,:,i);
    Z2s = Z2(:,:,i);
    
    %Z1o = Z1s;
    %Z1s = imadjust(Z1s);
    %Z2s = imadjust(Z2s);
    %     Zrgb = zeros(n,m,3);
    %     Zrgb(:,:,1) = Z1s;
    %     Zrgb(:,:,3) = Z2s;
    
    %Make image using color map
    Zrgb = MergeImageSpc(Z1s,Z2s,CMAP);
    Zrgb = flipud(Zrgb);
    %     [A,cm] = rgb2ind(Zrgb,256);
    
    %Additional Features
    %Zrgb = insertText(Zrgb, [10,10], 'Slavko is cool', 'FontSize',2, 'BoxColor', 'white' );
    %Zrgb = imresize(Zrgb, 2);
    
    % Generate plot with colormap and device body
    figure(1); clf;
    imt2 = image(X(1,:),flipud(Y(:,1)),Zrgb);
    hold on
    plotFrameGroup(frmgrp);
    set(gca,'ydir','normal')
    
    

    % Write PNGs and GIFs
    drawnow();
    if save_plots
        
        im = frame2im(getframe(gcf()));
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind,cm,[dirpath,savedir, '/PNG/',savedir,'_', num2str(i), '_PNG.png'])
        if i == 1
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif', 'Loopcount',Inf,'DelayTime',1);
        elseif i == size(Zstore,3)
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif','WriteMode','append','DelayTime',1);
        else
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif','WriteMode','append','DelayTime',1/frame_rate);
        end
        
    end
end