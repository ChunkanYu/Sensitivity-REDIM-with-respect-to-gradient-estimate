clear all;clc;
npsi=3;
% user specified system geometry
prompt = {'Enter the spatial domain L:',...
          'Enter the number of grid n_g:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','100'};
opts.Interpreter = 'tex';
system_parameter = inputdlg(prompt,dlgtitle,dims,definput,opts);
L =str2num(system_parameter{1}); 
ng=[str2num(system_parameter{2})];


% user specified system parameter
prompt = {'Enter the diffusion coefficient d:',...
          'Enter the rate coefficient k_1:',...
          'Enter the rate coefficient k_2:',...
          'Enter the rate coefficient k_3:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','1','10','20'};
opts.Interpreter = 'tex';
system_parameter = inputdlg(prompt,dlgtitle,dims,definput,opts);
D =str2num(system_parameter{1}); 
k=[str2num(system_parameter{2}),str2num(system_parameter{3}),...
    str2num(system_parameter{4})];


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

steady_state=y(end,:);

plot(x,steady_state(1,1:ng));hold on;
plot(x,steady_state(1,ng+1:2*ng));
plot(x,steady_state(1,2*ng+1:3*ng));
legend('\psi_1','\psi_2','\psi_3');
xlabel('{\it x}'); ylabel('\psi_i');

