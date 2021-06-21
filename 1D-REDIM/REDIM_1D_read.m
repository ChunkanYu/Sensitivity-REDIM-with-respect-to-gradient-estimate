function [REDIM_1D_info]=REDIM_1D_read(nspe,REDIM_1D_xxx,k)

% REDIM_1D_xxx=importdata(redim_string);
REDIM_1D_xxx=abs(REDIM_1D_xxx);
ng_redim=size(REDIM_1D_xxx,2)/(nspe+1);
gtheta=REDIM_1D_xxx(1,1:ng_redim);
REDIM_1D_state_xxx=REDIM_1D_xxx(1,ng_redim+1:end);

for i=1:nspe
    REDIM_1D_state(i,:)=REDIM_1D_state_xxx(1,(i-1)*ng_redim+1: i*ng_redim);
    dpsidtheta(i,:)=gradient(REDIM_1D_state(i,:),gtheta);
    d2psidtheta2(i,:)=gradient(dpsidtheta(i,:),gtheta);
end
REDIM_1D_info.state=REDIM_1D_state;
REDIM_1D_info.gtheta=gtheta;
REDIM_1D_info.dpsidtheta=dpsidtheta;
REDIM_1D_info.d2psidtheta2=d2psidtheta2;

for i=1:ng_redim
    dpsidtheta_perp(:,:,i)=gram_schmidt(dpsidtheta(:,i));
end
REDIM_1D_info.dpsidtheta_perp=dpsidtheta_perp;

for i=1:ng_redim
      F=Freac(REDIM_1D_state(:,i)',k)';
      G(:,i) = F;
end

REDIM_1D_info.G=G;


dGdpsi=[-k(1) 0 0;k(1) -k(2) 0; 0 k(2) -k(3)];

REDIM_1D_info.dGdpsi=dGdpsi;


end
