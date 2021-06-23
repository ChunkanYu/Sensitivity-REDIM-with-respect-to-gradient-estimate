clear all;clc;

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


prompt = {'Enter the number of grid n_g^1 in \theta_1 direction:',...
          'Enter the number of grid n_g^2 in \theta_2 direction:',...
          'Enter the gradient estimate \chi(\psi_1):'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'20','20','1'};
opts.Interpreter = 'tex';
system_parameter = inputdlg(prompt,dlgtitle,dims,definput,opts);
ng1 =str2num(system_parameter{1}); 
ng2 =str2num(system_parameter{2}); 
chi0=str2num(system_parameter{3});

fprintf('now the 2D REDIM will be generated, please wait ... \n');
[REDIM_2D]=REDIM_2D_generation(D,k,ng1,ng2,chi0);
fprintf('!!! now the 2D REDIM has been successfully generated!!! \n');

surf(REDIM_2D(:,:,3),REDIM_2D(:,:,4),REDIM_2D(:,:,5));

[sensitivity_psi2]=REDIM_2D_sensitivity(REDIM_2D,D,k,chi0);