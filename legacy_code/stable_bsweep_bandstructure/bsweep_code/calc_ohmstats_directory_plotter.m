for k=1:Nogs
    figure(k); hold on
    set(gca,'fontsize',14);
    set(gcf,'color','white');
    if sum(frm1.terminalOhmics==k)
        plot(B_exp,ohmNums(:,k));
        ylabel(['$I_{',num2str(k),'}/I_{in}$'],'interpreter','latex','FontSize',16,'fontweight','bold','FontName','Arial');
    elseif k==frm1.injector_ohmic
        plot(B_exp,V_source);
        ylabel('$V_{source}\ (arb.\ un.)$','interpreter','latex','FontSize',16,'FontName','Arial');
    else
        plot(B_exp,ohmFlux(:,k)./V_source);
        ylabel(['$V_{',num2str(k),'}/V_{source}$'],'interpreter','latex','FontSize',16,'FontName','Arial');
    end
   if exist('B','var')
        xlabel('$B (T)$','interpreter','latex','FontSize',16,'fontweight','bold','FontName','Arial');
   end
   set(gca,'xlim',[min(B_exp),max(B_exp)]);
%    xlim([-0.05,0.05])
         
end
%%
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

% k = 2;
% figure(k); hold on
% set(gca,'fontsize',14);
% set(gcf,'color','white');
% 
% plot(B_exp,ohmNumsp1(:,k)/2,'r');
% plot(B_exp,ohmNumsm1(:,k)/2,'b');
% ylabel('$I_{S4}/I_{in}$','interpreter','latex','FontSize',16,'fontweight','bold','FontName','Arial');
% 
% if exist('B','var')
%     xlabel('$B (T)$','interpreter','latex','FontSize',16,'fontweight','bold','FontName','Arial');
% end
% set(gca,'xlim',[min(B_exp),max(B_exp)]);
% %    xlim([-0.05,0.05])
% 
% 
% k = 4;
% figure(k); hold on
% set(gca,'fontsize',14);
% set(gcf,'color','white');
% 
% plot(B_exp,ohmNumsp1(:,k)/2,'r');
% plot(B_exp,ohmNumsm1(:,k)/2,'b');
% ylabel('$I_{F4}/I_{in}$','interpreter','latex','FontSize',16,'fontweight','bold','FontName','Arial');
% 
% if exist('B','var')
%     xlabel('$B (T)$','interpreter','latex','FontSize',16,'fontweight','bold','FontName','Arial');
% end
% set(gca,'xlim',[min(B_exp),max(B_exp)]);
% %    xlim([-0.05,0.05])
% 

