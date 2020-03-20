addpath('..')
dirpath = '..\..\Data\stable_bsweep_bandstructure\delafossite\Sherlock\';
savedir = 'rotData';
% subdir = '\xip1Data\';
subdir = '\';
savepath = [dirpath,savedir,subdir];
dirlist = ls([dirpath,savedir,subdir,'phi*']);
% ohmstatsPath = [dirpath,savedir,subdir];
ohmstatsPath = [dirpath,savedir];
Ndir = size(dirlist,1);

edgenumStore = 0;
N_injectStore = 0;
rngseeds = zeros(1,Ndir);

rngseedstore = 1;

for iDir = 1:Ndir
    dir = strtrim(squeeze(dirlist(iDir,:)));
    display(dir)
    
    fname = dir;
    
    load([savepath,fname,'\frameData.mat']);
    
    frm1=frmgrp.frms{1};
    
    addpath([savepath,fname])
    filelist=ls([savepath,fname,'\*classes.mat']);

    Nfiles=size(filelist,1);% number of files to loop through
    
    if exist('B','var')
        e = 1.60217662E-19; % C
        hbar = 1.054571800E-34; % Js
        B_exp = B * hbar/e * 10^16; % T
    else
        B_exp=1:Nfiles;
    end
    
    %% Load in all edge nums
    edgenums=zeros(Nfiles,frm1.Nedges);
    for i=1:Nfiles
        
        load([savepath,fname,'/',squeeze(filelist(i,:))]);
        
        edgenums(i,:)=edgenum;
        
    end
    
    %% Update stored parameters
    edgenumStore = edgenumStore + edgenums;
    N_injectStore = N_injectStore + N_inject;
    
    
    %% If storing rng seeds
    if rngseedstore
       load([savepath,fname,'/rngData.mat'])
       rngseeds(iDir) = seed;
    end

    %%
    fprintf('%i/%i\n',iDir,Ndir);
end

% Ndirp = Ndir;
% %%
% addpath('..')
% dirpath = '..\..\Data\stable_bsweep_bandstructure\miniBand\';
% % savedir = 'New folder';
% subdir = '\xim1Data\';
% savepath = [dirpath,savedir,subdir];
% dirlist = ls([dirpath,savedir,subdir,'phi*']);
% Ndir = size(dirlist,1);
% 
% rngseeds = zeros(1,Ndir);
% 
% for iDir = 1:min(Ndir,Ndirp)
%     dir = strtrim(squeeze(dirlist(iDir,:)));
%     display(dir)
%     
%     fname = dir;
%     
%     load([savepath,fname,'\frameData.mat']);
%     
%     frm1=frmgrp.frms{1};
%     
%     addpath([savepath,fname])
%     filelist=ls([savepath,fname,'\*classes.mat']);
% 
%     Nfiles=size(filelist,1);% number of files to loop through
%     
%     if exist('B','var')
%         e = 1.60217662E-19; % C
%         hbar = 1.054571800E-34; % Js
%         B_exp = B * hbar/e * 10^16; % T
%     else
%         B_exp=1:Nfiles;
%     end
%     
%     %% Load in all edge nums
%     edgenums=zeros(Nfiles,frm1.Nedges);
%     for i=1:Nfiles
%         
%         load([savepath,fname,'/',squeeze(filelist(i,:))]);
%         
%         edgenums(i,:)=edgenum;
%         
%     end
%     
%     %% Update stored parameters
%     edgenumStore = edgenumStore + edgenums;
%     N_injectStore = N_injectStore + N_inject;
%     
%     
%     %% If storing rng seeds
%     if rngseedstore
%        load([savepath,fname,'/rngData.mat'])
%        rngseeds(iDir) = seed;
%     end
% 
%     %%
%     fprintf('%i/%i\n',iDir,Ndir);
% end

%%
edgenums = edgenumStore;
N_inject = N_injectStore;


Nogs=length(frm1.ogs);
edgenumstot = zeros(Nfiles,Nogs);
ohmNums=zeros(Nfiles,Nogs);
ohmFlux=zeros(Nfiles,Nogs);

for k=1:Nogs
    edgenumstot(:,k)=sum(edgenums(:,frm1.ogs{k}.nedge_in_ohmic),2);
    ohmNums(:,k)=sum(edgenums(:,frm1.ogs{k}.nedge_in_ohmic),2)/N_inject;
    ohmFlux(:,k)=ohmNums(:,k)/frm1.ogs{k}.L;
end
V_source=ohmFlux(:,frm1.injector_ohmic) + 1/frm1.ogs{frm1.injector_ohmic}.L;

var = movvar(edgenumstot,11,1);

figure(1); clf; hold on;
xlabel('B (T)')
ylabel('Total collisions')
figure(2); clf; hold on;
xlabel('B (T)')
ylabel('N_{inject}*(Total Collisions)/(Source Collisions)')
for i = 1:8
    figure(1)
    plot(B_exp,edgenumstot(:,i))
    figure(2)
    plot(B_exp,N_inject*edgenumstot(:,i)./edgenumstot(:,8))
end
figure(3); clf;
plot(B_exp,var.^(1/2))
xlabel('B (T)')
ylabel('\sigma')

%%
mkdir([ohmstatsPath,'/ohmstats']);
save([ohmstatsPath,'/ohmstats/',savedir,'_','ohmstats.mat'],'edgenums','ohmNums','ohmFlux','V_source','N_inject');
save([ohmstatsPath,'/ohmstats/',savedir,'_','workspace.mat']);

%%
for k=1:Nogs
    figure(k)
    set(gca,'fontsize',14);
    set(gcf,'color','white');
    if sum(frm1.terminalOhmics==k)
        plot(B_exp,ohmNums(:,k));
        ylabel(['$I_{',num2str(k),'}/I_{in}$'],'interpreter','latex','FontSize',16,'fontweight','bold');
    elseif k==frm1.injector_ohmic
        plot(B_exp,V_source);
        ylabel('$V_{source}\ (arb.\ un.)$','interpreter','latex','FontSize',16);
    else
        plot(B_exp,ohmFlux(:,k)./V_source);
        ylabel(['$V_{',num2str(k),'}/V_{source}$'],'interpreter','latex','FontSize',16);
    end
   if exist('B','var')
        xlabel('$B (T)$','interpreter','latex','FontSize',16,'fontweight','bold');
   end
   set(gca,'xlim',[min(B_exp),max(B_exp)]);
%    xlim([-0.05,0.05])
   
   saveas(gcf(),[ohmstatsPath,'/ohmstats/',savedir,'_','ohmic_',num2str(k,'%02d'),'.fig'])
     im = frame2im(getframe(gcf()));
     [imind, cm] = rgb2ind(im, 256);
     imwrite(imind,cm,[ohmstatsPath,'/ohmstats/',savedir,'_','ohmic_',num2str(k,'%02d'),'.png'])  
       
end


figure(k+1)
frm1.plotOhmicLabels();
hold on;
plotFrameGroup(frmgrp);
hold off;
saveas(gcf(),[ohmstatsPath,'/ohmstats/',savedir,'_','ohmicsKey.fig'])
im = frame2im(getframe(gcf()));
[imind, cm] = rgb2ind(im, 256);
imwrite(imind,cm,[ohmstatsPath,'/ohmstats/',savedir,'_','ohmicsKey.png'])

%%
for k = 3:7
    figure(k+10)
    set(gca,'fontsize',14);
    set(gcf,'color','white');

        plot(B_exp,(ohmFlux(:,k)-ohmFlux(:,1))./V_source);
        ylabel(['$(V_{',num2str(k),'}-V_1)/V_{source}$'],'interpreter','latex','FontSize',16);
        xlabel('$B (T)$','interpreter','latex','FontSize',16,'fontweight','bold');

   set(gca,'xlim',[min(B_exp),max(B_exp)]);
%    xlim([-0.05,0.05])
   set(gca,'ylim',[-0.2,0.6]);
   
%    saveas(gcf(),[savepath,fname,'/ohmstats/',fname,'_','ohmic_',num2str(k,'%02d'),'.fig'])
%      im = frame2im(getframe(gcf()));
%      [imind, cm] = rgb2ind(im, 256);
%      imwrite(imind,cm,[savepath,fname,'/ohmstats/',fname,'_','ohmic_',num2str(k,'%02d'),'.png'])  
end

%%
for k=1:Nogs
    figure(k+20)
    set(gca,'fontsize',14);
    set(gcf,'color','white');
    if sum(frm1.terminalOhmics==k)
        plot(B_exp,ohmNums(:,k));
        ylabel(['$I_{',num2str(k),'}/I_{in}$'],'interpreter','latex','FontSize',16,'fontweight','bold');
    elseif k==frm1.injector_ohmic
        plot(B_exp,V_source);
        ylabel('$V_{source}\ (arb.\ un.)$','interpreter','latex','FontSize',16);
    else
        plot(B_exp,ohmFlux(:,k));
        ylabel(['$V_{',num2str(k),'}$'],'interpreter','latex','FontSize',16);
    end
    if exist('B','var')
        xlabel('$B (T)$','interpreter','latex','FontSize',16,'fontweight','bold');
    end
    set(gca,'xlim',[min(B_exp),max(B_exp)]);
    %     xlim([-0.05,0.05])
    
    saveas(gcf(),[ohmstatsPath,'/ohmstats/',savedir,'_','ohmic_',num2str(k,'%02d'),'.fig'])
    im = frame2im(getframe(gcf()));
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind,cm,[ohmstatsPath,'/ohmstats/',savedir,'_','ohmic_',num2str(k,'%02d'),'.png'])
end
