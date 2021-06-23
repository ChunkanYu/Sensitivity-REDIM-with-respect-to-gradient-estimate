function[yp]=REDIM_vector_field(t,y,odepar)
k=odepar.k; npsi=odepar.npsi; ng1=odepar.ng1; ng2=odepar.ng2;
D=odepar.D; gtheta1=odepar.gtheta1; gtheta2=odepar.gtheta2; 
chi0=odepar.chi0;

C=zeros(2,npsi); C(1,1)=1; C(2,3)=1;

for i=1:npsi
   xxx=y((i-1)*ng1*ng2+1:i*ng1*ng2); 
   for j=1:ng1
      psi(:,j,i)=xxx((j-1)*ng2+1:j*ng2);
   end
end

I=eye(npsi);
   
for i=1:npsi
    for j=1:ng2
       dpsidtheta1_1_matrix(j,:,i) = gradient(psi(j,:,i),gtheta1(j,:));
       d2psidtheta1_2_matrix(j,:,i) = gradient(dpsidtheta1_1_matrix(j,:,i),gtheta1(j,:));
    end
    for j=1:ng1
       dpsidtheta2_1_matrix(:,j,i) = gradient(psi(:,j,i),gtheta2(:,j));
       d2psidtheta2_2_matrix(:,j,i) = gradient(dpsidtheta2_1_matrix(:,j,i),gtheta2(:,j));
    end
end

chi_theta1=zeros(ng2,ng1); chi_theta2=zeros(ng2,ng1);
for i=2:ng1-1
    for j=2:ng2-1
       psitheta_now=[dpsidtheta1_1_matrix(j,i,1) dpsidtheta2_1_matrix(j,i,1);
                     dpsidtheta1_1_matrix(j,i,2) dpsidtheta2_1_matrix(j,i,2);
                     dpsidtheta1_1_matrix(j,i,3) dpsidtheta2_1_matrix(j,i,3)];
       
       chi_now=inv(C*psitheta_now)*C * [chi0 0 0]';
%        chi_now=pinv(psitheta_now) * [chi0 0 0]';
       chi_theta1(j,i)=chi_now(1,1);
       chi_theta2(j,i)=chi_now(2,1);
    end
end

for i=1:npsi
    for j=1:ng1
       d2psidtheta12_matrix(:,j,i) = gradient(dpsidtheta1_1_matrix(:,j,i),gtheta2(:,j));
    end
end



FF=zeros(ng2,ng1);
% chichi=[ng1 .01];
for i=2:ng1-1
    for j=2:ng2-1
        chichi=[chi_theta1(j,i) chi_theta2(j,i)];
        for kk=1:npsi
           dPsidTheta(kk,1:2)=[dpsidtheta1_1_matrix(j,i,kk),dpsidtheta2_1_matrix(j,i,kk)];
           d2PsidTheta2=[d2psidtheta1_2_matrix(j,i,kk),d2psidtheta12_matrix(j,i,kk);
                                d2psidtheta12_matrix(j,i,kk),d2psidtheta2_2_matrix(j,i,kk)];
           Diffusion(kk,1)=D*chichi*d2PsidTheta2*chichi';
           PSI_matrix(kk)=psi(j,i,kk);
        end
    invdPsidTheta=inv(C*dPsidTheta)*C;
    Proj=I-dPsidTheta*invdPsidTheta;
    
    F=Freac(PSI_matrix,k)';
    
    
    PHI=F+Diffusion;
    XX=Proj*PHI;
    for kk=1:npsi
        FF(j,i,kk)=XX(kk,1);
    end 
    end
end



yp=[];

for i=1:npsi
   FFF=FF(:,:,i);

   for j=1:ng1
       yp=[yp;FFF(:,j)];
   end
   
end

% sum(yp)

% surf(psi(:,:,1),psi(:,:,2),psi(:,:,3)); pause(.001);

   % plot(theta,y(ng_redim+1:2*ng_redim));hold on;pause(.001);

end

