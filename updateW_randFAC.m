function W=updateW_randFAC(X,Y,V,W0,beta,C)
[n,M]=size(Y);   
[J,Q]=size(V);
Y2=cen(Y);
numW=zeros(J,1);  
for i=1:n
    for j=1:J
        B(i,:,j)=V(j,:)'*X(i,j); 
    end
     Bi= reshape(B(i,:,:),Q,J);
    Expf=V'*W0*X(i,:)'+beta*((Y2(i,:))'-C'*V'*W0*X(i,:)');
    numW=numW+(Bi'*Expf);
end

denW=zeros(J,J);
for i=1:n
     for j=1:J
        B(i,:,j)=V(j,:)'*X(i,j); 
    end
     Bi= reshape(B(i,:,:),Q,J);
    denW=denW+Bi'*Bi;
end

W=diag(pinv(denW)*numW);


 
 

    