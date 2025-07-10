function  [Results]=CluRegApplication(rangeQ,nrandstart,X,Y) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% Code for application of real data set
% rangeQ is a vector. Example range: Q=[2,3,4,5]
% nrandstart number of starting points
% X matrix nxJ of covariates
% Y matrix nxM of responses
% V0 Variable-cluster membership matrix Initialization with kmeans 
% Output:
% Results list containing, for each value of Q, the best solution among the nrandstart initializations, corresponding to the following elements: 
% cpt Computation time 
% X matrix nxJ of covariates
% Y matrix nxM of responses
% Copt regression coefficient matrix
% Wopt Diagonal weight matrix 
% Vopt Variable-cluster membership matrix
% Sigmaeopt Covariance matrix
% LQmax Log-likelihood
% iteropt iterazion
% bicopt BIC value 
% aicopt AIC value 
% npar number of estimated parameters
% nrandstart number of starting points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=zscorepop(X);
Y=zscorepop(Y);
[n,M]=size(Y);
J=size(X,2);
LQmax=-inf; 
g=2;
Results={'rangeQ','cpt','X','Y','Copt','Wopt','Vopt','Sigmaeopt','LQmax','iteropt','BICopt','AICopt','npar','nrandstart'};

for i=1:size(rangeQ,2)
     disp(sprintf('Q=%g', rangeQ(i)))
      Q=rangeQ(i);
      t1=cputime;
for r=1:nrandstart
     V0=zeros(J,Q);  % Initialization kmeans 
      id=kmeans(X',Q,'Replicates',50,'EmptyAction','singleton');
      for j=1:J
          V0(j,id(j))=1;
      end 
     [V,W,C,SIGMAe,LQ,Ldif,iter,bic,aic,nparameters]=CluReg_randFAC(X,Y,Q,V0); % main function
      if Ldif<0
          Ldif, r=r-1;
      end
          
    if LQ>LQmax
        LQmax=LQ; Vopt=V; Wopt=W; Copt=C; 
        Sigmaeopt=SIGMAe; iteropt=iter; bicopt=bic; aicopt=aic; npar=nparameters;
    end   
end
 cpt=cputime-t1;
 disp(sprintf('cputime=%g', cpt));
Results{g,1}=rangeQ(i);
Results{g,2}=cpt;
Results{g,3}=X;
Results{g,4}=Y;
Results{g,5}=Copt;
Results{g,6}=Wopt;
Results{g,7}=Vopt;
Results{g,8}=Sigmaeopt;
Results{g,9}=LQmax;
Results{g,10}=iteropt;
Results{g,11}=bicopt;
Results{g,12}=aicopt;
Results{g,13}=npar;
Results{g,14}=nrandstart;
g=g+1;
end;


