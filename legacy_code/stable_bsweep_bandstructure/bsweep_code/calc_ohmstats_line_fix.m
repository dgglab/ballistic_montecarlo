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

%00
W = 6.7;
Y = 2*W-174/2;

%10
%W = 7.2;
%Y = 2*W-173/2;

%20
%W = 6.8;
%Y = 2*W-146/2;

%30
%W = 7.6;
%Y = 2*W-171/2;

xLine = [-4, 4];
yLine = [Y, Y]+1E-8;

for iDir = 1:Ndir
    tic 
    curDir = dirlist(iDir).name;
    display(curDir)
    
    
    frameData = load([savepath,curDir,'/frameData.mat']);
    frm1=frameData.frmgrp.frms{1};
    
    addpath([savepath,curDir])
    filelist = dir([savepath,curDir,'/*arcdata.txt']);
    
    classlist = dir([savepath,curDir,'/*classes.mat']);
    masterclasslist = dir([savepath,'masterclass','/*classes.mat']);
    Nfiles = length(filelist); % number of files to loop through
    
    linelist = dir([savepath,curDir,'/*lineohmstatsFIX.mat']);
    Nline = length(linelist);
    %     figure(n); clf; hold on;
    %     n=n+1;
    %     plot(frm1.px,frm1.py,'k','linewidth',1);
    
%     edgenumsLine=zeros(Nfiles,1);
    if Nline ~= Nfiles
        fprintf('running this directory\n')
        parfor iFile = 1:Nfiles 
            curFile = filelist(iFile).name;
            lineCheck = dir([savepath,curDir,'/',curFile,'_lineohmstatsFIX.mat']);
            check = length(lineCheck);
            if ~check
                display(curFile)
                fid=fopen([savepath,curDir,'/',curFile]);
                arcdata=fread(fid,'double');
                arcdata=reshape(arcdata,[],4);
                fclose(fid);

		% Generate save file
		edgenumsLine = nan;
		parsave_ohmline([savepath,curDir,'/',curFile,'_lineohmstatsFIX.mat'],edgenumsLine);        	

%         edgenums(i,:)=edgenum;

                % Load in acrclass
                curClass = classlist(iFile).name;
		mastercurClass = masterclasslist(iFile).name;
                classes = load([savepath,'masterclass','/',mastercurClass]);
        
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
            %         drawnow
                end
%         edgenumsLine(iFile) = edgenumStore;
                edgenumsLine = edgenumStore;
                fprintf('ran suc:%s\n',curFile);
                parsave_ohmline([savepath,curDir,'/',curFile,'_lineohmstatsFIX.mat'],edgenumsLine);
            end
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

