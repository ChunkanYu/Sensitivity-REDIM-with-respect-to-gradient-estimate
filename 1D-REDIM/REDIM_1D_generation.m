function [REDIM_1D]=REDIM_1D_generation(D,k,ng,chi0)
npsi=3;


SS=importdata('SS.mat');

ns_SS=size(SS,2)/npsi;

inipro(1,:)=linspace(SS(1,1),SS(1,ns_SS),ng);

for i=2:3
   inipro(1,(i-1)*ng+1:i*ng) = interp1(SS(1,1:ns_SS),SS((i-1)*ns_SS+1:i*ns_SS),inipro(1,1:ng)); 
end


gtheta=1:1:ng;

odepar.k=k; odepar.npsi=npsi; odepar.ng=ng;
odepar.D=D;  odepar.chi0=chi0;
odepar.gtheta= gtheta;
opts = odeset('reltol',1e-4,'abstol',1e-6); 

[t,y] = ode15s(@REDIM_vector_field,[0,100],inipro',opts,odepar); % 15s

REDIM_1D=[gtheta y(end,:)];
