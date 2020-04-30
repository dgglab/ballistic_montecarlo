addpath('..')

dirpath = '../../Data/stable_bsweep_bandstructure/';
savedir = 'delafossiteBar';
% subdir = '\xip1Data\';
subdir = '/';
savepath = [dirpath,savedir,subdir];
dirlist = dir([dirpath,savedir,subdir,'phi*'])
% ohmstatsPath = [dirpath,savedir,subdir];
ohmstatsPath = [dirpath,savedir];
Ndir = length(dirlist);


dirpath = '..\..\Data\stable_bsweep_bandstructure\delafossiteBar\';
savedir = 'phi_pi6';
% subdir = '\xip1Data\';
subdir = '\';
savepath = [dirpath,savedir,subdir];
dirlist = ls([dirpath,savedir,subdir,'phi*']);
% ohmstatsPath = [dirpath,savedir,subdir];
ohmstatsPath = [dirpath,savedir];
Ndir = size(dirlist,1);


W = 6.7;
Y = 2*W-174/2;

xLine = [-4, 4];
yLine = [Y, Y];
tic
pause(1800*rand)
toc

for iDir = 1:Ndir
    tic 
    curDir = dirlist(iDir).name;
    display(curDir)
    
    
    frameData = load([savepath,curDir,'/frameData.mat']);
    frm1=frameData.frmgrp.frms{1};
    
    addpath([savepath,curDir])
    filelist = dir([savepath,curDir,'/*arcdata.txt']);
    
    classlist = dir([savepath,curDir,'/*classes.mat']);
    Nfiles = length(filelist); % number of files to loop through
    
    linelist = dir([savepath,curDir,'/*lineohmstats.mat']);
    Nline = length(linelist);


for iDir = 1:Ndir
    tic
    
    dir = strtrim(squeeze(dirlist(iDir,:)));
    display(dir)
    fname = dir;
    
    frameData = load([savepath,fname,'\frameData.mat']);
    frm1=frameData.frmgrp.frms{1};
    
    addpath([savepath,fname])
    filelist = ls([savepath,fname,'/*arcdata.txt']);
    
    classlist = ls([savepath,fname,'/*classes.mat']);
    Nfiles=size(filelist,1);% number of files to loop through
    
    
    % Check if we have already run this directory
    linelist = ls([savepath,fname,'/*lineohmstats.mat']);
    

    %     figure(n); clf; hold on;
    %     n=n+1;
    %     plot(frm1.px,frm1.py,'k','linewidth',1);
    

%     edgenumsLine=zeros(Nfiles,1);
    if Nline ~= Nfiles
    %if Nline == 0
        fprintf('running this directory\n')
        parfor iFile = 1:Nfiles 
            curFile = filelist(iFile).name;
            fid=fopen([savepath,curDir,'/',curFile]);
            arcdata=fread(fid,'double');
            arcdata=reshape(arcdata,[],4);
            fclose(fid);
        
%         edgenums(i,:)=edgenum;

            % Load in acrclass
            curClass = classlist(iFile).name;
            classes = load([savepath,curDir,'/',curClass]);
        
            Narcs=size(arcdata,1);
    
            edgenumStore = 0;
%         flag = 0;
            for iArc=1:Narcs
%             arcdata(iArc,:)
%             x1 = arcdata(1);
%             y1 = arcdata(2);
%             ind1 = arcdata(3);
%             ind2 = arcdata(4);
                [xout,yout] = classes.arcclass.arc(arcdata(iArc,:));
                [xi,yi] = polyxpoly(xLine,yLine,xout,yout);
             

    %     edgenumsLine=zeros(Nfiles,1);
    if strcmp(linelist,'')
        parfor iFile=1:Nfiles
            
            fid=fopen([savepath,fname,'/',squeeze(filelist(iFile,:))]);
            arcdata=fread(fid,'double');
            arcdata=reshape(arcdata,[],4);
            fclose(fid);
            
            %         edgenums(i,:)=edgenum;
            
            % Load in acrclass
            classes = load([savepath,fname,'/',squeeze(classlist(iFile,:))]);
            
            Narcs=size(arcdata,1);
            
            edgenumStore = 0;
            %         flag = 0;
            %for iArc=1:Narcs
            for iArc=1:2
                %             arcdata(iArc,:)
                %             x1 = arcdata(1);
                %             y1 = arcdata(2);:wq
                %             ind1 = arcdata(3);
                %             ind2 = arcdata(4);
                [xout,yout] = classes.arcclass.arc(arcdata(iArc,:));
                [xi,yi] = polyxpoly(xLine,yLine,xout,yout);
                

                if mod(length(xi),2) == 0
                    edgenumTemp = 0;
                else
                    edgenumTemp = 1;
                end

            
%             edgenumTemp = length(xi);
%             if ~isempty(xi)
%                 edgenumTemp = 1;
%                 if flag == 1
%                     edgenumTemp = 0;
%                 end
%                 flag = 1;
%             else
%                 edgenumTemp = 0;
%                 flag = 0;
%             end
                edgenumStore = edgenumStore + edgenumTemp;
            
            %         in = inpolygon(xout,yout,frm1.px,frm1.py);
            %         if size(xout(~in),2)>2
            %             badarcnum(badcount) = i;
            %             badarcs(badcount,:) = arcdata(i,:);
            %             badclass{badcount} = arcclass;
            %             badcount = badcount + 1;
            %         end
            %         plot(xout,yout,'.-')
                drawnow
            end
%         edgenumsLine(iFile) = edgenumStore;
            edgenumsLine = edgenumStore;

            parsave_ohmline([savepath,curDir,'/',curFile,'_lineohmstats.mat'],edgenumsLine);
        end
    %     parsave([savepath,fname,'/','lineohmstats.mat']);

                
                %             edgenumTemp = length(xi);
                %             if ~isempty(xi)
                %                 edgenumTemp = 1;
                %                 if flag == 1
                %                     edgenumTemp = 0;
                %                 end
                %                 flag = 1;
                %             else
                %                 edgenumTemp = 0;
                %                 flag = 0;
                %             end
                edgenumStore = edgenumStore + edgenumTemp;
                
                %         in = inpolygon(xout,yout,frm1.px,frm1.py);
                %         if size(xout(~in),2)>2
                %             badarcnum(badcount) = i;
                %             badarcs(badcount,:) = arcdata(i,:);
                %             badclass{badcount} = arcclass;
                %             badcount = badcount + 1;
                %         end
                %         plot(xout,yout,'.-')
                drawnow
            end
            %         edgenumsLine(iFile) = edgenumStore;
            edgenumsLine = edgenumStore;
            
            parsave_ohmline([savepath,fname,'/',squeeze(filelist(iFile,:)),'_lineohmstats.mat'],edgenumsLine);
        end
        %     parsave([savepath,fname,'/','lineohmstats.mat']);

        fprintf('%i/%i\n',iDir,Ndir);
        fprintf('finished running directory\n')
        toc
        break
    else
        fprintf('directory already contains linestats\n')
    end
end


