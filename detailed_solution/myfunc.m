function yp=myfunc(t,PSI,odepar)
M=odepar.M; npsi=odepar.npsi; k=odepar.k; ng=odepar.ng;
PSI_matrix=zeros(ng,npsi);


for i=1:npsi
    PSI_matrix(:,i)=PSI((i-1)*ng+1: i*ng);
end
for i=1:size(PSI_matrix,1)
    F(i,:)=Freac(PSI_matrix(i,:),k);
end

for i=1:size(PSI_matrix,2)
    FF(:,i) = F(:,i) + M*PSI_matrix(:,i);
end

FF(1,:)=0; FF(end,:)=0;
yp=[];
for i=1:npsi
    yp=[yp;FF(:,i)];
end

end

function F=Freac(psi,k)
psi1=psi(1); psi2=psi(2); psi3=psi(3);
k1=k(1); k2=k(2); k3=k(3);
F(1,1)=-k1*psi1;
F(1,2)=k1*psi1-k2*psi2;
F(1,3)=k2*psi2-k3*psi3;
end
