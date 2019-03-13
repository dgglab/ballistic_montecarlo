addpath('..')
dirpath = '..\..\Data\stable_bsweep_bandstructure\delafossiteBar\';
savedir = 'bar30_run7';
% subdir = '\xip1Data\';
subdir = '\';
savepath = [dirpath,savedir,subdir];
dirlist = ls([dirpath,savedir,subdir,'phi*']);
% ohmstatsPath = [dirpath,savedir,subdir];
ohmstatsPath = [dirpath,savedir];
Ndir = size(dirlist,1);


for iDir = 1:Ndir
    tic
    
    dir = strtrim(squeeze(dirlist(iDir,:)));
    display(dir)
    fname = dir;
    
    frameData = load([savepath,fname,'\frameData.mat']);
    frm1=frameData.frmgrp.frms{1};
    
    addpath([savepath,fname])
%     filelist = ls([savepath,fname,'/*arcdata.txt']);
    filelist = ls([savepath,fname,'/*lineohmstatsFIX.mat']);
    
    classlist = ls([savepath,'masterclass','/*classes.mat']);
    
    Nfiles=size(filelist,1);% number of files to loop through
    edgenumStore = zeros(Nfiles,1);
    
    %     figure(n); clf; hold on;
    %     n=n+1;
    %     plot(frm1.px,frm1.py,'k','linewidth',1);
    
%     edgenumsLine=zeros(Nfiles,1);
    for iFile=1:Nfiles
        
        ohmStats = load([savepath,fname,'/',squeeze(filelist(iFile,:))]);
        edgenumStore(iFile) = ohmStats.edgenumLine;


        
    end
%     parsave([savepath,fname,'/','lineohmstats.mat']);
    fprintf('%i/%i\n',iDir,Ndir);
    
    toc

end

edgenumStore = (edgenumStore-N_inject)./2 + N_inject;