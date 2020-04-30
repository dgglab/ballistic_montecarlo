function [outx, outy]=rotate_and_offset(inx,iny,theta,delx,dely)

outx=cos(theta)*inx+sin(theta)*iny+delx;
outy=cos(theta)*iny-sin(theta)*inx+dely;