function [REDIM_2D_info]=REDIM_2D_read(nspe,REDIM_2D_xxx,k)

ng1=size(REDIM_2D_xxx,2); ng2=size(REDIM_2D_xxx,1);

REDIM_2D_xxx=abs(REDIM_2D_xxx);
gtheta1=REDIM_2D_xxx(:,:,1); gtheta2=REDIM_2D_xxx(:,:,2);
REDIM_2D_state_xxx=REDIM_2D_xxx(:,:,3:end);

REDIM_2D_info.state=REDIM_2D_state_xxx;
REDIM_2D_info.gtheta1=gtheta1;
REDIM_2D_info.gtheta2=gtheta2;

for kk=1:nspe
    %
    for j=1:ng2
      dpsidtheta1_1_matrix(j,:,kk)=gradient(REDIM_2D_state_xxx(j,:,kk),gtheta1(j,:));
      d2psidtheta1_2_matrix(j,:,kk)=gradient(dpsidtheta1_1_matrix(j,:,kk),gtheta1(j,:));
    end
    %
    for i=1:ng1
      dpsidtheta2_1_matrix(:,i,kk)=gradient(REDIM_2D_state_xxx(:,i,kk),gtheta2(:,i));
      d2psidtheta2_2_matrix(:,i,kk)=gradient(dpsidtheta2_1_matrix(:,i,kk),gtheta2(:,i));
    end
    %
    for i=1:ng1
        d2psidtheta12_matrix(:,i,kk)=gradient(dpsidtheta1_1_matrix(:,i,kk),gtheta2(:,i));
    end
    %
end

REDIM_2D_info.dpsidtheta1_1_matrix=dpsidtheta1_1_matrix;
REDIM_2D_info.dpsidtheta2_1_matrix=dpsidtheta2_1_matrix;
REDIM_2D_info.d2psidtheta1_2_matrix=d2psidtheta1_2_matrix;
REDIM_2D_info.d2psidtheta2_2_matrix=d2psidtheta2_2_matrix;
REDIM_2D_info.d2psidtheta12_matrix=d2psidtheta12_matrix;


for i=1:ng1
   for j=1:ng2
      a=dpsidtheta1_1_matrix(j,i,1);b=dpsidtheta1_1_matrix(j,i,2);c=dpsidtheta1_1_matrix(j,i,3);
      d=dpsidtheta2_1_matrix(j,i,1);f=dpsidtheta2_1_matrix(j,i,2);g=dpsidtheta2_1_matrix(j,i,3);
      vvv(1,1)=(b*g-c*f)/(a*f-b*d);
      vvv(2,1)=(a*g-c*d)/(b*d-a*f);
      vvv(3,1)=1;
      dpsidtheta_perp(:,j,i)=vvv/norm(vvv);
   end
end

REDIM_2D_info.dpsidtheta_perp=dpsidtheta_perp;

for i=1:ng1
   for j=1:ng2
      state_now=[REDIM_2D_state_xxx(j,i,1);REDIM_2D_state_xxx(j,i,2);REDIM_2D_state_xxx(j,i,3)];
      F=Freac(state_now',k)';
      G(:,j,i) = F;
   end
end

REDIM_2D_info.G=G;


dGdpsi=[-k(1) 0 0;k(1) -k(2) 0; 0 k(2) -k(3)];

REDIM_2D_info.dGdpsi=dGdpsi;


end