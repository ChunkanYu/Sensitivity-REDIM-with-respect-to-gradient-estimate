function[yp]=REDIM_vector_field(t,y,odepar)
npsi=odepar.npsi; D=odepar.D; k=odepar.k;
ng=odepar.ng; gtheta=odepar.gtheta; chi0=odepar.chi0;
ng_redim=ng;
I=eye(npsi);
PSI=y;
   
   PSI_matrix=zeros(ng,npsi);
   for i=1:npsi
      PSI_matrix(:,i)=PSI((i-1)*ng_redim+1: i*ng_redim);
   end
   
for i=1:npsi
    dpsidtheta_matrix(:,i) = gradient(PSI_matrix(:,i),gtheta);
    d2psidtheta2_matrix(:,i) = gradient(dpsidtheta_matrix(:,i),gtheta);
end


C=zeros(1,npsi); C(1,1)=1;


for i=1:ng
    dpsidtheta_now=[dpsidtheta_matrix(i,1);
                    dpsidtheta_matrix(i,2);
                    dpsidtheta_matrix(i,3)];
   chi(i)=inv(C*dpsidtheta_now)*C*[chi0;0;0]; 
end



for i=1:ng
    for j=1:npsi
        dPsidTheta(j,1) = interp1(gtheta,dpsidtheta_matrix(:,j),gtheta(i));
        d2PsidTheta2(j,1) = interp1(gtheta,d2psidtheta2_matrix(:,j),gtheta(i));
        psi(j,1) = interp1(gtheta,PSI_matrix(:,j),gtheta(i));
    end
    
    invdPsidTheta=inv(C*dPsidTheta)*C;
    Proj=I-dPsidTheta*invdPsidTheta;
    
    F=Freac(PSI_matrix(i,:),k)';
    PHI=F+D*chi(1,i)*chi(1,i)*d2PsidTheta2; 
%     size(Proj)
%     size(PHI)
    FF(:,i)=Proj*PHI;
end

FF(:,1)=FF(:,1)*0; FF(:,ng_redim)=FF(:,ng_redim)*0;

FF=FF';

yp=[];
   for i=1:npsi
       yp=[yp;FF(:,i)];
   end

end


