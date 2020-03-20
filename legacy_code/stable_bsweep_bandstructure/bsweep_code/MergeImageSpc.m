function [ IM ] = MergeImageSpc( I1,I2,CMAP)
%Author - Slavko Rebec (email: srebec@stanford.edu)
%Description - Creates an image using two intensity maps and a premade
%              color map from the function Colormap2d

%Inputs:
%I1 - Intensity map with range inside of [0,1], of same size of I2
%I2 - Intensity map with range inside of [0,1], of same size of I1
%CMAP - Color map made using Colormap2d

%Output:
%IM - RGB image of size I1 x 3 

%Remaps range of I1 and I2 to CMAP indexes
[x1,x2, ~] = size(CMAP);
I1 = round(I1*(x1-1)+1);
I2 = round(I2*(x2-1)+1);

%Initializes Output Image
[n,m] = size(I1);
IM = zeros(n,m,3);

%Indexes entries from CMAP using I1 and I2
a = I1(:);
b = I2(:);
c = sub2ind([x1,x2],a,b);
CMAP_1 = CMAP(:,:,1);
CMAP_2 = CMAP(:,:,2);
CMAP_3 = CMAP(:,:,3);
s1 = CMAP_1(c);
s1 = reshape(s1,n,m);
s2 = CMAP_2(c);
s2 = reshape(s2,n,m);
s3 = CMAP_3(c);
s3 = reshape(s3,n,m);
IM(:,:,1) = s1;
IM(:,:,2) = s2;
IM(:,:,3) = s3;

%Shorter but slower
% inda = sub2ind([x1,x2],I1(:),I2(:));
% IMa = [CMAP(inda);CMAP(inda+x1*x2);CMAP(inda+2*x1*x2)];
% IM = reshape(IMa, n,m,3);

end
