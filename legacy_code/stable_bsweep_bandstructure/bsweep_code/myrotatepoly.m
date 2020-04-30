function [outx, outy]=myrotatepoly(inx,iny,theta)

outx=cos(theta)*inx+sin(theta)*iny;
outy=cos(theta)*iny-sin(theta)*inx;