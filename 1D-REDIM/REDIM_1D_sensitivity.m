% clear all;clc;
function [sensitivity_vector,Term,MM,yend]=REDIM_1D_sensitivity(REDIM_1D,D,k,chi0)
ptheta=1; % ng=150;  
npsi=3;

I=eye(npsi); C=I(:,ptheta); C_perp=[I(:,ptheta+1:end)];

% REDIM_1D_xxx=importdata('inipro.mat');

[REDIM_1D_info]=REDIM_1D_read(npsi,REDIM_1D,k);

ng_redim= size(REDIM_1D_info.gtheta,2);

chi=chi0*gradient(REDIM_1D_info.gtheta,REDIM_1D_info.state(1,:));


I=eye(ng_redim);
% A(theta)
for i=1:ng_redim
    A_xx=inv(REDIM_1D_info.dpsidtheta_perp(:,:,i)'*C_perp)*REDIM_1D_info.dpsidtheta_perp(:,:,i)';
    A_theta(:,:,i)=A_xx*REDIM_1D_info.dGdpsi*C_perp;
    clear A_xx;
end

% B(theta)
for i=1:ng_redim
    B_xx = inv(C'*REDIM_1D_info.dpsidtheta(:,i))*C'; 
    B_theta(i,1)=B_xx * (REDIM_1D_info.G(:,i)+D*chi(i)*chi(i)*REDIM_1D_info.d2psidtheta2(:,i));
    clear B_xx;
end

% C(theta)
for i=1:ng_redim
    C_xx=inv(REDIM_1D_info.dpsidtheta_perp(:,:,i)'*C_perp);
    C_xxx=C_xx*REDIM_1D_info.dpsidtheta_perp(:,:,i)'*REDIM_1D_info.d2psidtheta2(:,i);
    C_theta(:,i)=2*D*C_xxx*chi(i)*chi(i);
    clear C_xx; clear C_xxx;
end


Matrix_II=zeros(ng_redim,ng_redim);
Matrix_III=Matrix_II;
Matrix_IV=zeros(ng_redim,npsi-1);
Matrix_a11=Matrix_II;Matrix_a12=Matrix_II;Matrix_a21=Matrix_II;Matrix_a22=Matrix_II;
Matrix_II(1,1)=1; Matrix_II(ng_redim,ng_redim)=1;
for i=2:ng_redim-1
    Matrix_II(i,i-1) = - 0.5*B_theta(i,1)/1;
    Matrix_II(i,i+1) = - Matrix_II(i,i-1);
    Matrix_III(i,i-1) = D * chi(i)*chi(i)/1;
    Matrix_III(i,i) = -2*Matrix_III(i,i-1);
    Matrix_III(i,i+1) = Matrix_III(i,i-1);
    for j=1:npsi-1
       Matrix_IV(i,j)=-C_theta(j,i); 
    end
    Matrix_a11(i,i)=A_theta(1,1,i);
    Matrix_a12(i,i)=A_theta(1,2,i);
    Matrix_a21(i,i)=A_theta(2,1,i);
    Matrix_a22(i,i)=A_theta(2,2,i);
end

Matrix_I=[Matrix_a11 Matrix_a12;    Matrix_a21 Matrix_a22];
Matrix_II_now=[Matrix_II Matrix_II*0; 0*Matrix_II Matrix_II];
Matrix_III_now=[Matrix_III Matrix_III*0; 0*Matrix_III Matrix_III];
Matrix_IV_now=Matrix_IV(:,1);
for j=2:npsi-1
   Matrix_IV_now=[Matrix_IV_now;Matrix_IV(:,j)]; 
end

MM=Matrix_I+Matrix_II_now+Matrix_III_now;

yendend=MM\Matrix_IV_now;
yend=[yendend(1:ng_redim,1)';yendend(ng_redim+1:2*ng_redim,1)'];

for i=1:npsi-1
   dsigmadtheta(i,:)=gradient(yend(i,:),REDIM_1D_info.gtheta);
   d2sigmadtheta2(i,:)=gradient(dsigmadtheta(i,:),REDIM_1D_info.gtheta);
end


for i=1:ng_redim
   sensitivity(:,i)=C_perp * yend(:,i); 
   Term_1(:,i) = A_theta(:,:,i)*yend(:,i);
   Term_2(:,i) = dsigmadtheta(:,i)*B_theta(i,1);
   Term_3(:,i) = D*chi(i)*chi(i)*d2sigmadtheta2(:,i);
   Term_4(:,i) = C_theta(:,i);
end


sensitivity_vector=reshape(sensitivity',[],1)';

Term.Term_1 = Term_1;Term.Term_2 = Term_2;
Term.Term_3 = Term_3;Term.Term_4 = Term_4;

end
