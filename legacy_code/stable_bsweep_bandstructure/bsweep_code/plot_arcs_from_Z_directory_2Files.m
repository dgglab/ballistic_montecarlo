% Set to true to export png files
save_plots=true;
% Frame rate of gif output
frame_rate=5;

Contrast=5;
skip = 1;

fixPlot = 1; % Manually removes points outside the boundary

%%
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

Zstore1 = Zstore;

%%

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

cmap1=ones(128,3);
cmap1(:,1)=linspace(1,0.85,128);
cmap1(:,2)=linspace(1,0,128);
cmap1(:,3)=linspace(1,0,128);

cmap2=ones(128,3);
cmap2(:,1)=linspace(1,0,128);
cmap2(:,2)=linspace(1,0.4470,128);
cmap2(:,3)=linspace(1,0.7410,128);

colormap([cmap1;cmap2])
climmax = (N_inject/Contrast)*255;



for i = 1:size(Zstore,3)
    
    figure(1); clf;
    set(gcf,'color','w');
    set(gcf,'color','white');
    imt1 = image(X(1,:),Y(:,1),min(Zstore1(:,:,i),N_inject/Contrast)/(N_inject/Contrast)*255);
    set(gca,'ydir','normal')
    hold on
    plotFrameGroup(frmgrp);
    colormap(cmap.red);
    
    figure(2); clf;
    set(gcf,'color','white');
    imt2 = image(X(1,:),Y(:,1),min(Zstore2(:,:,i),N_inject/Contrast)/(N_inject/Contrast)*255);
    set(gca,'ydir','normal')
    hold on
    plotFrameGroup(frmgrp);
    colormap(cmap.blue);
    
    
    imtest1 = frame2im(getframe(figure(1)));
    imtest2 = frame2im(getframe(figure(2)));

        imadded = uint8(min(double(imtest1),double(imtest2)));
    figure(5); clf; hold on
    set(gcf,'color','white');
    imshow(imadded)
% %     
%     figure(3); clf; hold on
%      set(gcf,'color','white');
%     imf = imfuse(imtest1,imtest2,'ColorChannels',[2 1 1]);
%     ima = imshow(imf);


    
%     mesh = 255*ones(size(Zstore1(:,:,i)));
%     figure(1); clf; hold on
%     im1new = image(X(1,:),Y(:,1),cat(3,mesh,0*mesh,0*mesh));
%     set(gca,'ydir','normal')
%     
%     
%     im2new = image(X(1,:),Y(:,1),cat(3,0*mesh,0*mesh,1*mesh));
%     set(gca,'ydir','normal')
%     plotFrameGroup(frmgrp);
% %     Zalpha1 = Zstore1(:,:,i)/max(max(Zstore1(:,:,i)));
% %     Zalpha2 = Zstore2(:,:,i)/max(max(Zstore2(:,:,i)));
%     Zalpha1 = min(Zstore1(:,:,i),N_inject/Contrast)/(N_inject/Contrast);
%     Zalpha2 = min(Zstore2(:,:,i),N_inject/Contrast)/(N_inject/Contrast);
% 
% 
%     set(im1new,'AlphaData',0.5*(Zalpha1))
%     set(im2new,'AlphaData',0.5*(Zalpha2))

    drawnow();
    if save_plots
        
        im = frame2im(getframe(gcf()));
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind,cm,[dirpath,savedir, '/PNG/',savedir,'_', num2str(i), '_PNG.png'])
        if i == 1
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif', 'Loopcount',Inf,'DelayTime',1/frame_rate);
        else
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif','WriteMode','append','DelayTime',1/frame_rate);
        end
        
    end
end