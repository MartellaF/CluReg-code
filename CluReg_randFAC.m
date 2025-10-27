%%%%%%%%%%%%%%%%%%%%%% Main function %%%%%%%%%%%%%%%%%%%%%%%%%%
% Multivariate Regression Model Based on Latent
% Predictors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,W,C,SIGMAe,LQ,Ldif,iter,bic,aic,nparameters]=CluReg_randFAC(X,Y,Q,V0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X (n x J) Explanatory variables
% V (J x Q) membership matrix for clustering variables
% Y (n x M) Response variables
% V (J x Q) binary and row stochastic
% W (J x J) diagonal weight matrix with (W*V) column orthogonal
% C (Q x M) regression coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=size(V0,1); 
[n,M]=size(Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization phase 
eps     = 0.1;  
maxiter = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Starting points                                                      
W0=diag(rand(J,1)); 
% Normalization
Vtild0=W0*V0; 
Vtild0=Vtild0*pinv(Vtild0'*Vtild0)^0.5;
W0=diag(sum(Vtild0*pinv(V0),2)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y2=cen(Y);
F0=X*W0*V0;                                                 
C0=pinv(F0'*F0)*F0'*(Y2);                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIGMAe0=diag(diag((cov(Y2-F0*C0))));           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta0=C0*inv(C0'*C0+SIGMAe0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L0=-inf;
%ltol=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration phase
for iter=1:maxiter
%%%%%%%%%%%%%%%%%%%%%%%%%%% Update C 
numC=zeros(Q,M);
for i=1:n
    Expf=V0'*W0*X(i,:)'+beta0*((Y2(i,:)')-C0'*V0'*W0*X(i,:)') ;
  numC=numC+Expf*(Y2(i,:));
end
denC=zeros(Q,Q);
for i=1:n
    Expff=(eye(Q)-beta0*C0')+((V0'*W0*X(i,:)'+beta0*((Y2(i,:))'-C0'*V0'*W0*X(i,:)'))*(V0'*W0*X(i,:)'+beta0*((Y2(i,:))'-C0'*V0'*W0*X(i,:)'))');
    denC=denC+Expff;
end
C=inv(denC)*numC;                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%% Update W
W=updateW_randFAC(X,Y,V0,W0,beta0,C0);  
%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization W
  Vtild=W*V0; 
  Vtild=Vtild*pinv(Vtild'*Vtild)^0.5;
  W=diag(sum(Vtild*pinv(V0),2));
%%%%%%%%%%%%%%%%%%%%%%%%%%% Update V  
  V=zeros(J,Q);
  IDE=eye(Q);
  Vapp=V0; 
      for j=1:J
         tt=zeros(Q,1);
         for v=1:Q
            V0(j,:)=IDE(v,:); 
            if all(sum(V0))==0, tt(v)=-inf; break, end
               tt(v)=Likeli1_randFAC(X,Y,SIGMAe0,W,V0,C); 
         end
         [~,pos]=max(tt); 
         V(j,pos)=1;
         V0(j,:)=V(j,:);
            if sum((sum(V0)==0))>=1  %'a cluster with 0 observation'
                 V(j,:)=Vapp(j,: );
                  V0(j,:)=Vapp(j,: );
           end
      end    
secondsum=zeros(M,Q);
for i=1:n
    Expf=V0'*W0*X(i,:)'+beta0*((Y2(i,:))'-C0'*V0'*W0*X(i,:)') ;
  secondsum= secondsum+Y2(i,:)'*Expf';
end;     
%%%%%%%%%%%%%%%%%%%%%%%%%%% Update SIGMAe 
SIGMAe=1/n*diag(diag(((Y2)'*(Y2)-secondsum*C)));
%%%%%%%%%%%%%%%%%%%%%%%%%%% Update beta
beta=C*pinv(C'*C+SIGMAe);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute likelihood
[LQ]=Likeli1_randFAC(X,Y,SIGMAe,W,V,C);
Ldif = LQ-L0;
 if   Ldif > eps 
     L0=LQ;
     V0=V;
     C0=C;
     W0=W;
     SIGMAe0=SIGMAe;
     beta0=beta;
     if (iter>maxiter) 
            break
     end
 elseif Ldif < 0
              break
 elseif               Ldif < eps    
           break   
 end;
end 
nparameters=J+(J-Q)+(M*Q)+M;
bic=-2*LQ+log(n)*nparameters;
aic=-LQ+nparameters;



