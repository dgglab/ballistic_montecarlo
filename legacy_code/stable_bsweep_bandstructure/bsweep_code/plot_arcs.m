
dirpath = '..\..\Data\stable_bsweep_bandstructure\delafossite\Sherlock\';
savedir = 'rotData\';
savepath = [dirpath,savedir];
dirlist = ls([dirpath,savedir,'\phi*']);
addpath([savepath,dir])
load([savepath,dir,'\frameData.mat']);

%find the files associated with simulation results
fname = dir;
filelist = ls([savepath,fname,'/*arcdata.txt']);

classlist = ls([savepath,fname,'/*classes.mat']);

Nfiles = size(filelist,1);% number of files to loop through



%generate X,Y meshgrid for binning arc points
% save([savepath,fname,'/imageProperties.mat'],'X','Y','dL');


n=1;
badcount=1;
%for k = 1:Nfiles
for k = 868
    tic
    figure(n); hold on;
    n=n+1;
    %plot(frm1.px,frm1.py,'k','linewidth',1);
    plotFrameGroup(frmgrp);
    hold on
    % Load in data file
    fid=fopen([savepath,fname,'/',squeeze(filelist(k,:))]);
    arcdata=fread(fid,'double');
    arcdata=reshape(arcdata,[],4);
    fclose(fid);
    
    % Load in acrclass
    load([savepath,fname,'/',squeeze(classlist(k,:))]);
      
    Narcs=size(arcdata,1);
    
    
    %for i=1650:1655
    %for i=1259
%     for i = [768:772, 1259, 1650:1655]
    %for i = 1650:Narcs
     for i = [768:772, 1259, 11856:11861]
        arcdata(i,:)
        [xout,yout] = arcclass.arc(arcdata(i,:));
        plotFrameGroup(frmgrp);
        hold on
        plot(xout,yout,'k-','Linewidth',2)
        drawnow
        
%        if xout(1) > 1.712 && xout(1) < 2.538
%            if yout(1) <-5.15 && yout(1) > -5.35
%                i
%                figure(1); clf; hold on;
%                plotFrameGroup(frmgrp);
%                hold on
%                for j = 0:20
%                    [xout,yout] = arcclass.arc(arcdata(i+j,:));
%                    plot(xout,yout,'.-')
%                end
%                drawnow
%            end
%        end
        
        %in = inpolygon(xout,yout,frm1.px,frm1.py);
        %if size(xout(~in),2)>2
        %    badarcnum(badcount) = i;
        %    badarcs(badcount,:) = arcdata(i,:);
        %    badclass{badcount} = arcclass;
        %    badcount = badcount + 1;
        %end

    end
    drawnow
    %save(['./',fname,'/',fname,'_',num2str(k,'%03d'),'Z.mat'],'Z');
    k
    toc
    

    
    


end
%end