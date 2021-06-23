% clear all;clc;
function [sensitivity_psi2]=REDIM_2D_sensitivity(REDIM_2D,D,k,chi0)
%  REDIM_2D=importdata('REDIM_2D_xxx.mat'); D=1; k=[ 1 10 20]; chi0=1;
npsi=3;

I=eye(npsi); C=I(:,[1 3]); C_perp=I(:,2);


[REDIM_2D_info]=REDIM_2D_read(npsi,REDIM_2D,k);

ng1=size(REDIM_2D,2); ng2=size(REDIM_2D,1);

chi_gtheta1=zeros(ng2,ng1); chi_gtheta2=zeros(ng2,ng1);
for i=1:ng1
    for j=1:ng2
       psitheta_now=[REDIM_2D_info.dpsidtheta1_1_matrix(j,i,1) REDIM_2D_info.dpsidtheta2_1_matrix(j,i,1);
                     REDIM_2D_info.dpsidtheta1_1_matrix(j,i,2) REDIM_2D_info.dpsidtheta2_1_matrix(j,i,2);
                     REDIM_2D_info.dpsidtheta1_1_matrix(j,i,3) REDIM_2D_info.dpsidtheta2_1_matrix(j,i,3)];
       
       chi_now=inv(C'*psitheta_now)*C' * [chi0 0 0]';
       chi_gtheta1(j,i)=chi_now(1,1);
       chi_gtheta2(j,i)=chi_now(2,1);
    end
end

% A(theta)
for j=1:ng2
    for i=1:ng1
    A_xx=inv(REDIM_2D_info.dpsidtheta_perp(:,j,i)'*C_perp)*REDIM_2D_info.dpsidtheta_perp(:,j,i)';
    A_theta(j,i)=A_xx*REDIM_2D_info.dGdpsi*C_perp;
    clear A_xx;
    end
end

% B(theta)
for j=1:ng2
    for i=1:ng1
        chichi=[chi_gtheta1(j,i) chi_gtheta2(j,i)];
        psitheta=[REDIM_2D_info.dpsidtheta1_1_matrix(j,i,1),REDIM_2D_info.dpsidtheta2_1_matrix(j,i,1);
                  REDIM_2D_info.dpsidtheta1_1_matrix(j,i,2),REDIM_2D_info.dpsidtheta2_1_matrix(j,i,2);
                  REDIM_2D_info.dpsidtheta1_1_matrix(j,i,3),REDIM_2D_info.dpsidtheta2_1_matrix(j,i,3);];
     B_xx = inv(C'*psitheta)*C'; 
     for kk=1:npsi
            psithetatheta=[REDIM_2D_info.d2psidtheta1_2_matrix(j,i,kk), REDIM_2D_info.d2psidtheta12_matrix(j,i,kk);
                           REDIM_2D_info.d2psidtheta12_matrix(j,i,kk), REDIM_2D_info.d2psidtheta2_2_matrix(j,i,kk)];
            B_xxx=chichi*psithetatheta*chichi'; B_xxxx(kk,1)=D*B_xxx;
            clear C_xxx;
     end
     B_theta(:,j,i)= -B_xx*(REDIM_2D_info.G(:,j,i) + B_xxxx);
     clear B_xx;
    end
end

% C(theta)
for j=1:ng2
    for i=1:ng1
        chichi=[chi_gtheta1(j,i) chi_gtheta2(j,i)];
        C_xx=inv(REDIM_2D_info.dpsidtheta_perp(:,j,i)'*C_perp)*REDIM_2D_info.dpsidtheta_perp(:,j,i)';
        for kk=1:npsi
            psithetatheta=[REDIM_2D_info.d2psidtheta1_2_matrix(j,i,kk), REDIM_2D_info.d2psidtheta12_matrix(j,i,kk);
                           REDIM_2D_info.d2psidtheta12_matrix(j,i,kk), REDIM_2D_info.d2psidtheta2_2_matrix(j,i,kk)];
            C_xxx(kk,1)=chichi*psithetatheta*chichi';
        end
        C_theta(j,i) = 2*D*C_xx*C_xxx;
        clear C_xx; clear C_xxx;
    end
end

sigma=zeros(ng2,ng1);
sigma0=[];
for i=1:ng1
    sigma0=[sigma0;sigma(:,i)];
end

opts = odeset('reltol',1e-4,'abstol',1e-6); 

odepar.ng1=ng1; odepar.ng2=ng2;odepar.D=D; 
odepar.chi_gtheta1=chi_gtheta1; odepar.chi_gtheta2=chi_gtheta2;
odepar.A_theta=A_theta; odepar.B_theta=B_theta;odepar.C_theta=C_theta;
odepar.gtheta1=REDIM_2D_info.gtheta1;
odepar.gtheta2=REDIM_2D_info.gtheta2;

[t,y] = ode15s(@REDIM_2D_sensitivity_func,[0,1],sigma0,opts,odepar); % 15s

yendend=y(end,:);

for i=1:ng1
      yend(:,i)=yendend((i-1)*ng2+1:i*ng2); 
   end

for j=1:ng2
    for i=1:ng1
       aaa=C_perp * yend(j,i);
       sensitivity_psi2(j,i)=aaa(2,1);
    end
end

surf(REDIM_2D_info.gtheta1,REDIM_2D_info.gtheta2,sensitivity_psi2);

end