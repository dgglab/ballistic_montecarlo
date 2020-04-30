
dirpath = '..\..\Data\stable_bsweep_bandstructure\delafossiteBar\';
savedir = 'bar30_run3\';
savepath = [dirpath,savedir];
dirlist = ls([dirpath,savedir,'\phi*']);
addpath([savepath,dir])
load([savepath,dir,'\frameData.mat']);

%find the files associated with simulation results
fname = dir;
filelist = ls([savepath,fname,'/*arcdata.txt']);

classlist = ls([savepath,'masterclass','/*classes.mat']);

Nfiles = size(filelist,1);% number of files to loop through



%generate X,Y meshgrid for binning arc points
% save([savepath,fname,'/imageProperties.mat'],'X','Y','dL');

filenum = 310;

filelist = ls([savepath,fname,'/*',num2str(filenum),'arcdata.txt']);

classlist = ls([savepath,'masterclass','/*',num2str(filenum),'classes.mat']);

Nfiles = size(filelist,1);% number of files to loop through


W = 6.7;
Y = 2*W-174/2;

xLine = [-4, 4];
yLine = [Y, Y]+1E-8;

n=11;
badcount=1;
for k = 1:Nfiles
    tic
    figure(n); clf; hold on;
    n=n+1;
    plot(frm1.px,frm1.py,'k','linewidth',1);
    line(xLine,yLine)
    % Load in data file
    fid=fopen([savepath,fname,'/',squeeze(filelist(k,:))]);
    arcdata=fread(fid,'double');
    arcdata=reshape(arcdata,[],4);
    fclose(fid);
    
    % Load in acrclass
    load([savepath,'masterclass','/',squeeze(classlist(k,:))]);
      
    Narcs=size(arcdata,1);
    
    
    for i=1:Narcs
        %i
        %arcdata(i,:)
        [xout,yout] = arcclass.arc(arcdata(i,:));
        
%         in = inpolygon(xout,yout,frm1.px,frm1.py);
%         if size(xout(~in),2)>2
%             badarcnum(badcount) = i;
%             badarcs(badcount,:) = arcdata(i,:);
%             badclass{badcount} = arcclass;
%             badcount = badcount + 1;
%         end

        [xi,yi] = polyxpoly(xLine,yLine,xout,yout);

        if length(xi)>0
            disp(length(xi))
        end

        if mod(i,1000)==0
            figure(n-1); clf; hold on;
            plot(frm1.px,frm1.py,'k','linewidth',1);
            line(xLine,yLine)
            plot(xout,yout,'.-')
            drawnow
        else
            plot(xout,yout,'.-')
            drawnow
        end
    end
    drawnow
    %save(['./',fname,'/',fname,'_',num2str(k,'%03d'),'Z.mat'],'Z');
    k
    toc
    

    
    


end
%end