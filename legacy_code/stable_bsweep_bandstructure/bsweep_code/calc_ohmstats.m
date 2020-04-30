load([savepath,fname,'/frameData.mat']);

frm1=frmgrp.frms{1};

addpath([savepath,fname])
filelist=ls([savepath,fname,'/*classes.mat']);
k=1;
while k<=size(filelist,1)
    if  ~isempty(strfind(filelist(k,:),'/*Z.mat'))
        filelist(k,:)=[];
    else
        k=k+1;
    end
end
Nfiles=size(filelist,1);% number of files to loop through

if exist('B','var')
    e = 1.60217662E-19; % C
    hbar = 1.054571800E-34; % Js
    B_exp = B * hbar/e * 10^16; % T
else
    B_exp=1:Nfiles;
end

%%
edgenums=zeros(Nfiles,frm1.Nedges);
for i=1:Nfiles
    
    load([savepath,fname,'/',squeeze(filelist(i,:))]);
   
    edgenums(i,:)=edgenum;

end

% Nogs=length(frm1.ogs);
Nogs = 2;
ohmNums=zeros(Nfiles,Nogs);
ohmFlux=zeros(Nfiles,Nogs);

for k=1:Nogs
    ohmNums(:,k)=sum(edgenums(:,frm1.ogs{k}.nedge_in_ohmic),2)/N_inject;
    ohmFlux(:,k)=ohmNums(:,k)/frm1.ogs{k}.L;
end
V_source=ohmFlux(:,frm1.injector_ohmic) + 1/frm1.ogs{frm1.injector_ohmic}.L;

%%
mkdir([savepath,fname,'/ohmstats']);

save([savepath,fname,'/ohmstats/',fname,'_','ohmstats.mat'],'edgenums','ohmNums','ohmFlux','V_source');
%%



for k=1:Nogs
    figure(k)
    set(gca,'fontsize',14);
    set(gcf,'color','white');
    if sum(frm1.terminalOhmics==k);
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
   
   saveas(gcf(),[savepath,fname,'/ohmstats/',fname,'_','ohmic_',num2str(k,'%02d'),'.fig'])
     im = frame2im(getframe(gcf()));
     [imind, cm] = rgb2ind(im, 256);
     imwrite(imind,cm,[savepath,fname,'/ohmstats/',fname,'_','ohmic_',num2str(k,'%02d'),'.png'])  
       
end


figure(k+1)
frm1.plotOhmicLabels();
hold on;
plotFrameGroup(frmgrp);
hold off;
saveas(gcf(),[savepath,fname,'/ohmstats/',fname,'_','ohmicsKey.fig'])
im = frame2im(getframe(gcf()));
[imind, cm] = rgb2ind(im, 256);
imwrite(imind,cm,[savepath,fname,'/ohmstats/',fname,'_','ohmicsKey.png'])


   
