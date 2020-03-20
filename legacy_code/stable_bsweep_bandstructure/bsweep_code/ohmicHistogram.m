addpath('..')

dirpath = '../../Data/stable_bsweep_bandstructure/';
savedir = 'delafossiteBar/';
% subdir = '\xip1Data\';
subdir = 'bar30_run6/';
savepath = [dirpath,savedir,subdir];
dirlist = ls([dirpath,savedir,subdir,'phi*']);
% ohmstatsPath = [dirpath,savedir,subdir];
ohmstatsPath = [dirpath,savedir];
Ndir = size(dirlist,1);


for iDir = 1:Ndir
    
    curDir = strtrim(squeeze(dirlist(iDir,:)));
    display(curDir)
    fname = curDir;
    
    frameData = load([savepath,fname,'\frameData.mat']);
    frm1=frameData.frmgrp.frms{1};
    
    addpath([savepath,fname])
    filelist = ls([savepath,fname,'/*arcdata.txt']);
    
    classlist = ls([savepath,fname,'/*classes.mat']);
    Nfiles=size(filelist,1);% number of files to loop through
    
    
    % Check if we have already run this directory
    linelist = ls([savepath,fname,'/*lineohmstats.mat']);
    injectInd = NaN(1E5,Nfiles);

        for iFile = 1:Nfiles 
            curFile = strtrim(squeeze(filelist(iFile,:)));
            display(curFile)
            fid=fopen([savepath,curDir,'/',curFile]);
            arcdata=fread(fid,'double');
            arcdata=reshape(arcdata,[],4);
            fclose(fid);
            injected = arcdata(:,3);
            injectInd(1:length(injected),iFile) = injected;
%         edgenums(i,:)=edgenum;

            % Load in acrclass
            curClass = strtrim(squeeze(classlist(iDir,:)));
            classes = load([savepath,curDir,'/',curClass]);
        
            Narcs=size(arcdata,1);
    
        end
end

%% 
for i = 1:5:size(injectInd,2)
    figure(i);clf;
    histogram(injectInd(:,i),253)
end