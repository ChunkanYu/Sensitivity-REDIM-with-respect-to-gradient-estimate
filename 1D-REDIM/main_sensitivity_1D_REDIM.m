clear all;clc; 
% This is the main code to generate the sensitivity of 1D REDIM with respect to gradient estimate.
% Several input options are needed, which will be required by GUI-Window
% The default values are consistent to the case studied in the paper

% user specified system parameter
prompt = {'Enter the diffusion coefficient d:',...
          'Enter the rate coefficient k_1:',...
          'Enter the rate coefficient k_2:',...
          'Enter the rate coefficient k_3:'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'1','1','10','20'};
opts.Interpreter = 'tex';
system_parameter = inputdlg(prompt,dlgtitle,dims,definput,opts);
D =str2num(system_parameter{1}); 
k=[str2num(system_parameter{2}),str2num(system_parameter{3}),...
    str2num(system_parameter{4})];


prompt = {'Enter the number of grid n_g:',...
          'Enter the gradient estimate \chi(\psi_1):'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'100','1'};
opts.Interpreter = 'tex';
system_parameter = inputdlg(prompt,dlgtitle,dims,definput,opts);
ng =str2num(system_parameter{1}); 
chi0=str2num(system_parameter{2});



[REDIM_1D]=REDIM_1D_generation(D,k,ng,chi0);

[sensitivity_vector,Term]=REDIM_1D_sensitivity(REDIM_1D,D,k,chi0);

state=REDIM_1D(ng+1:end);
state_psi1=state(1:ng);
figure(1);
plot(state_psi1,sensitivity_vector(ng+1:2*ng)); hold on;
plot(state_psi1,sensitivity_vector(2*ng+1:3*ng));

figure(2)
plot(state_psi1,Term.Term_1(1,:),'k-');hold on;
plot(state_psi1,Term.Term_2(1,:),'b-');
plot(state_psi1,Term.Term_3(1,:),'g-');
plot(state_psi1,Term.Term_4(1,:),'m-');

figure(3)
plot(state_psi1,Term.Term_1(2,:),'k-');hold on;
plot(state_psi1,Term.Term_2(2,:),'b-');
plot(state_psi1,Term.Term_3(2,:),'g-');
plot(state_psi1,Term.Term_4(2,:),'m-');

