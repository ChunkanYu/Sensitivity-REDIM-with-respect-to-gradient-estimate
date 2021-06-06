clear all;clc; close all;
L=5; D=1; npsi=3;
k=[1 10 20];
ng=100;
x=linspace(0,L,ng); dx=x(2)-x(1);
M=zeros(ng,ng);
for i=2:ng-1
    M(i,i-1) = 1; M(i,i)=-2; M(i,i+1)=1; end
M=M*D/dx/dx;

IC=zeros(ng,npsi); IC(1,1)=1;

ICC=reshape(IC,[],1);

odepar.M=M; odepar.k=k; odepar.npsi=npsi; odepar.ng=ng;


opts = odeset('reltol',1e-5,'abstol',1e-7); 

[t,y] = ode15s(@myfunc,[0,1e+2],ICC,opts,odepar); % 15s



for i=1:size(y,1)
   yyy=reshape(y(i,:),[],npsi);
   
    plot(x,yyy(:,1)); hold on;plot(x,yyy(:,2)*10);plot(x,yyy(:,3)*10);
    legend('\psi_!','10 \cdot \psi_2','10 \cdot \psi_3');
    xlabel('{\it t} /s'); ylabel('\psi');
   title(num2str(t(i,1))); pause(.1); hold off;
end