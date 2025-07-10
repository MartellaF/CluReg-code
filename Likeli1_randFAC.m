function [LQ]=Likeli1_randFAC(X,Y,SIGMAe,W,V,C)
  Y2=cen(Y);          
% 2 - 02 - 2018
[n,M]=size(Y); 
A=W*V*C; 
SIGMAy=(C'*C+SIGMAe);
SYinv=pinv(SIGMAy);  
deterY=det(SIGMAy);
costante=-0.5*n*M*log(2*pi);
Q=-0.5*trace(SYinv*(Y2-X*A)'*(Y2-X*A)); 
LQ=Q+costante-0.5*n*log(deterY);
