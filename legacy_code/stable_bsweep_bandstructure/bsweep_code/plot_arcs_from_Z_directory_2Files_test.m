%%

frame_rate=5;

Contrast=3;
skip = 1;

save_plots = 0;
if save_plots
    mkdir([dirpath,savedir, '/PNG'])
    mkdir([dirpath,savedir, '/GIF'])
end

cmap1=ones(128,3);
cmap1(:,1)=linspace(1,0.85,128);
cmap1(:,2)=linspace(1,0,128);
cmap1(:,3)=linspace(1,0,128);

cmap2=ones(128,3);
cmap2(:,1)=linspace(1,0,128);
cmap2(:,2)=linspace(1,0.4470,128);
cmap2(:,3)=linspace(1,0.7410,128);

colormap([cmap1;cmap2])
climmax = (N_inject/Contrast);

myCmaps;

%%
in = inpolygon(X,Y,frmgrp.px,frmgrp.py);
Zstore1fix = Zstore1;
Zstore2fix = Zstore2;
for i = 1:size(Zstore1,3)
    Zstore1fix(:,:,i) = Zstore1(:,:,i).*in;
    Zstore2fix(:,:,i) = Zstore2(:,:,i).*in;
end

%%

for i = 1:size(Zstore,3)
%     
%     figure(1); clf; hold on
%     set(gcf,'color','w');
%     Z1 = Zstore1(:,:,i);
%     Z1(Z1>climmax) = climmax;
%     Z1 = 1+127*Z1/climmax;
%     x = [minX,maxX];
%     y = [minY,maxY];
%     
%     %im1 = imagesc(fliplr(Y(:,1)),(X(1,:)),Z1);
%     im1 = imagesc(x,y,Z1);
%     
%     
% %     im1 = imagesc(X(:,1),Y(1,:),Z1);
%     shading flat;
%     %im1.AlphaData = .5;
% 
%     Z2 = Zstore2(:,:,i);
%     Z2(Z2>climmax) = climmax;
%     Z2 = 129+127*Z2/climmax;
%     im2 = imagesc(x,y,Z2);
%     shading flat;
%     im2.AlphaData = .5;
% 
%   %using image. Z is thresholded and normalized to 255
    figure(2)
    set(gcf,'color','none');
    imt1 = image(X(1,:),Y(:,1),min(Zstore1fix(:,:,i),N_inject/Contrast)/(N_inject/Contrast)*255);
    set(gca,'ydir','normal')
    hold on
    plotFrameGroup(frmgrp);
    %use cmap of your choice
    colormap(cmap.red);
    
    figure(3)
    set(gcf,'color','none');
    imt2 = image(X(1,:),Y(:,1),min(Zstore2fix(:,:,i),N_inject/Contrast)/(N_inject/Contrast)*255);
    set(gca,'ydir','normal')
    hold on
    
    %use cmap of your choice
    colormap(cmap.blue);
    plotFrameGroup(frmgrp);
    
    figure(4)
    plotFrameGroup(frmgrp);
    
    imtest1 = frame2im(getframe(figure(2)));
    imtest2 = frame2im(getframe(figure(3)));
    imtest3 = frame2im(getframe(figure(4)));
  %  imtest1 = ind2rgb(imt1,imt1.CDataMapping.red);
    %imtest2 = ind2rgb(imt2,cmap.blue);
    
    

    %imaddeli = imfuse(imtest1,imtest2,'Scaling','joint','ColorChannels',[2 1 1]);
%     imadded = uint8(mod(double(imtest1)+double(imtest2),256));
    %imadded = mod(double(imtest1) + double(imtest2),256);
%     imadded = cat(3,imtest1(:,:,1), imtest1(:,:,2), imtest2(:,:,3));
%     imadded = uint8(min(mod(double(imtest1)+double(imtest2),256),255));
%     imadded = imcomplement(imtest1 + imtest2 + imcomplement(imtest3));
    %imadded = imcomplement(imcomplement(imtest1) + imcomplement(imtest2));
    
    %imadded = imcomplement(uint8(double(imcomplement(imtest1)) + double(imcomplement(imtest2))));
    %imadded = (uint8(max(double((imtest1)), double((imtest2)))));
%     imadded = (uint8((double((imtest1))+double((imtest2))-255)));
    %imadded = uint8(double((imtest1)) + double((imtest2)));
    imadded = uint8(min(double(imtest1),double(imtest2)));
    figure(5); clf; hold on
    imshow(imadded)
    
%     imf = imfuse(imtest1,imtest2,'ColorChannels',[2 1 1]);
%     ima = imshow(imf);
%     alphamap = double(sum(imadded,3)) ~=0;
%     set(ima,'AlphaData',alphamap)
%     hold on
%     
    %plot the frame and ohmics for visual reference

    
    
    
%     mesh = 255*ones(size(Zstore1(:,:,i)));
%     figure(1); clf; hold on
%     im1new = image(X(1,:),Y(:,1),cat(3,mesh,0*mesh,0*mesh));
%     set(gca,'ydir','normal')
%     
%     
%     im2new = image(X(1,:),Y(:,1),cat(3,0*mesh,0*mesh,1*mesh));
%     set(gca,'ydir','normal')
%     
% %     Zalpha1 = Zstore1(:,:,i)/max(max(Zstore1(:,:,i)));
% %     Zalpha2 = Zstore2(:,:,i)/max(max(Zstore2(:,:,i)));
%     Zalpha1 = min(Zstore1(:,:,i),N_inject/Contrast)/(N_inject/Contrast);
%     Zalpha2 = min(Zstore2(:,:,i),N_inject/Contrast)/(N_inject/Contrast);
% 
% 
%     set(im1new,'AlphaData',0.5*(Zalpha1))
%     set(im2new,'AlphaData',0.5*(Zalpha2))
%     


% 
%     
%     %using image. Z is thresholded and normalized to 255
%     image(X(1,:),Y(:,1),min(Zstore(:,:,i),N_inject/Contrast)/(N_inject/Contrast)*255);
%     hold on
%     
    %plot the frame and ohmics for visual reference
    
%     %use cmap of your choice
%     colormap(cmap.red);
    drawnow();
    if save_plots
        
        im = frame2im(getframe(gcf()));
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind,cm,[dirpath,savedir, '/PNG/',savedir,'_', num2str(i), '_PNG.png'])
        if i == 1
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif', 'Loopcount',Inf,'DelayTime',1/frame_rate);
        else
            imwrite(imind,cm,[dirpath,savedir, '/GIF/',savedir,'.gif'],'gif','WriteMode','append','DelayTime',1/frame_rate);
        end
        
    end
end