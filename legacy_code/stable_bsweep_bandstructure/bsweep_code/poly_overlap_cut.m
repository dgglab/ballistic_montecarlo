function [px,py, edgestats1]=poly_overlap_cut(px1,py1,px2,py2,edgestats1,edgecode)

N1=length(px1);
N2=length(px2);


[xcut,ycut,iicut]=polyxpoly(px1,py1,px2,py2);


%sorts to enable iterative insertion
A=[xcut,ycut,iicut];
A=sortrows(A,3);

 i=1;
 
 while(i<size(A,1)); 
     Nrepeat=size((find(A(:,3)==A(i,3))),1);
     Atemp=A(i-1+(1:Nrepeat),:);
     Lsqr=zeros(Nrepeat,1);
     for j=1:Nrepeat
         Lsqr(j)=(Atemp(j,1)-px1(A(i,3))).^2+(Atemp(j,2)-py1(A(i,3))).^2;
     end
     Atemp=[Atemp,Lsqr];
     Atemp=sortrows(Atemp,5);
     A(i-1+(1:Nrepeat),:)=Atemp(1:Nrepeat,1:4);
     i=i+Nrepeat;
end

xcut=A(:,1);
ycut=A(:,2);
iicut=A(:,3:4);


    
Ncuts=size(iicut,1);
N_on_edge=zeros(1,N1+Ncuts);

for i=1:Ncuts
    px1=[px1(1:(iicut(i,1)+i-1)),xcut(i),px1((iicut(i,1)+i):(N1+i-1))];
    py1=[py1(1:(iicut(i,1)+i-1)),ycut(i),py1((iicut(i,1)+i):(N1+i-1))];
    %if(~mod(iicut(i,1)+i-1,2))
    edgestats1=[edgestats1(1:(iicut(i,1)+i-1)),edgestats1((iicut(i,1)+i-1)),edgestats1((iicut(i,1)+i):(N1+i-2))];
    %else
    %edgestats1=[edgestats1(1:(iicut(i,1)+i-2)),edgecode,edgestats1((iicut(i,1)+i-1):(N1+i-2))];
    %end    
    N_on_edge(iicut(i,1)+i)=1;
end
    
px=px1;
py=py1;

IN=inpolygon(px,py,px2,py2);
%includes all internal and edge points
IN=IN|N_on_edge;

a=find(IN);

internal_edges=a(find(1==a(2:length(a))-a(1:(length(a)-1))));
edgestats1(internal_edges)=edgecode;





