clf
for k=0:5:100
    fname=['e_h_test_',num2str(k,'%03d')];
for i=1:51
    
    
    tic;
    
    load(['./',fname,'/',fname,'_',num2str(i,'%03d'),'.mat']);
   % load(['bsweep_Vprobe_33percent_scatter_',num2str(i,'%03d'),'.mat'],'edgenum');toc
    edgenums(i,:)=edgenum;
    toc
end



%%

xscale=0.035;
yscale=1E4/3;  
subplot(3,1,1)
hold on
R_1 = (sum(edgenums(:,frm1.ogs{4}.nedge_in_ohmic),2) / frm1.ogs{4}.L -sum(edgenums(:,frm1.ogs{3}.nedge_in_ohmic),2) / frm1.ogs{3}.L)/yscale;
plot(R_1+k/10,'linewidth',2)
set(gca,'ylim',[-1 10])
set(gca,'fontsize',14,'fontweight','bold')
set(gca,'xtick',[])
ylabel('R_1 (arb. un.)')
subplot(3,1,2)
hold on
R_2 = (sum(edgenums(:,frm1.ogs{7}.nedge_in_ohmic),2) / frm1.ogs{7}.L-sum(edgenums(:,frm1.ogs{3}.nedge_in_ohmic),2) / frm1.ogs{3}.L)/yscale;
plot(R_2+k/10,'linewidth',2)
set(gca,'ylim',[-1 10])
set(gca,'fontsize',14,'fontweight','bold')
ylabel('R_2 (arb. un.)')
subplot(3,1,3)
hold on
R_3 = (sum(edgenums(:,frm1.ogs{1}.nedge_in_ohmic),2) / frm1.ogs{1}.L-sum(edgenums(:,frm1.ogs{7}.nedge_in_ohmic),2) / frm1.ogs{7}.L)/yscale;
plot(R_3+k/10,'linewidth',2)
set(gca,'ylim',[-1 15])
set(gca,'fontsize',14,'fontweight','bold')
set(gcf,'color','white')    

%pk_ht_ratio(k+1) = R_3(4)/R_3(8);
%pk_area_ratio(k+1) = trapz(R_3(3:5))/trapz(R_3(6:10));
end
%figure
%subplot(2,1,1)
%plot(0:0.1:1,pk_ht_ratio)
%ylabel('Peak Height Ratio (R_{pk1} / R_{pk2})')
%subplot(2,1,2)
%plot(0:0.1:1,pk_area_ratio)
%ylabel('Peak Area Ratio')
%xlabel('p_{scatter}')

