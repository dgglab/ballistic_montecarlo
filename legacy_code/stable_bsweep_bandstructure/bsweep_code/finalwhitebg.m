clear
close all
clc


load('matlab (1).mat')
xr = X(1,:);
xL = xr(end) - xr(1);
yr = Y(:,1);
yL = yr(end) - yr(1);
g = 1/3;
Z1m = max(Zstore1(:));
Z2m = max(Zstore1(:));
Zmax = max([Z1m,Z2m]);
[n,m,ims] = size(Zstore1);

rgb = [1,0,.8,.4;1,0,0,0;1,.8,0,.4];
%rgb = [1,0,.85,.42;1,.447,.325,0.37;1,.741,0,.37];
sfina  = CustomColormap(50, 50, rgb );

for i = 1:ims
    Z1 = Zstore1(:,:,i);
    Z2 = Zstore2(:,:,i);
%     Z1s = 0.15*log(Z1+1);
%     Z2s = 0.15*log(Z2+1);
    Z1s = (Z1./Zmax*1.1).^g;
    Z2s = (Z2./Zmax*1.1).^g;
    Z1o = Z1s;
    %Z1s = imadjust(Z1s);

    %Z2s = imadjust(Z2s);
%     Zrgb = zeros(n,m,3);
%     Zrgb(:,:,1) = Z1s;
%     Zrgb(:,:,3) = Z2s;
    Zrgb = MergeImageSpc(Z1s,Z2s,sfina);
    Zrgb = flipud(Zrgb);
    [A,cm] = rgb2ind(Zrgb,256);
    
    if i == 1
        imwrite(A,cm,'test.gif','gif','LoopCount',Inf,'DelayTime',2);
    elseif i == ims
        imwrite(A,cm,'test.gif','gif','WriteMode','append','DelayTime',2);
    else
        imwrite(A,cm,'test.gif','gif','WriteMode','append','DelayTime',0.1);
    end
end


