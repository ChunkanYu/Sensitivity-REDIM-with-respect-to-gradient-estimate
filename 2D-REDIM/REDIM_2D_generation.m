function [REDIM_2D]=REDIM_2D_generation(D,k,ng1,ng2,chi0)
npsi=3;

[y0,psi_all]=PREPREDIM_read(npsi,ng1,ng2);

gtheta1=1:1:ng1; gtheta1=repmat(gtheta1,ng2,1);
gtheta2=[1:1:ng2]'; gtheta2=repmat(gtheta2,1,ng1);


odepar.k=k; odepar.npsi=npsi; odepar.ng1=ng1; odepar.ng2=ng2;
odepar.D=D; odepar.gtheta1=gtheta1; odepar.gtheta2=gtheta2; 
odepar.chi0=chi0;

opts = odeset('reltol',1e-4,'abstol',1e-6); 

[t,y] = ode15s(@REDIM_vector_field,[0,1],y0,opts,odepar); % 15s

yend=reshape(y(end,:),[],npsi);
psi_1=reshape(yend(:,1),ng2,ng1);
psi_2=reshape(yend(:,2),ng2,ng1);
psi_3=reshape(yend(:,3),ng2,ng1);

REDIM_2D(:,:,1)= gtheta1; REDIM_2D(:,:,2)= gtheta2;
for i=1:npsi
REDIM_2D(:,:,i+2)= reshape(yend(:,i),ng2,ng1);
end


