function [dy]=REDIM_2D_sensitivity_func(t,y,odepar)
   ng1=odepar.ng1; ng2=odepar.ng2;D=odepar.D; 
chi_gtheta1=odepar.chi_gtheta1; chi_gtheta2=odepar.chi_gtheta2;
A_theta=odepar.A_theta; B_theta=odepar.B_theta;C_theta=odepar.C_theta;
gtheta1=odepar.gtheta1; gtheta2=odepar.gtheta2;

   for i=1:ng1
      sigma(:,i)=y((i-1)*ng2+1:i*ng2); 
   end
   
   for j=1:ng2
      dsigmadtheta1_1_matrix(j,:)=gradient(sigma(j,:),gtheta1(j,:));
      d2sigmadtheta1_2_matrix(j,:)=gradient(dsigmadtheta1_1_matrix(j,:),gtheta1(j,:));
    end
    %
    for i=1:ng1
      dsigmadtheta2_1_matrix(:,i)=gradient(sigma(:,i),gtheta2(:,i));
      d2sigmadtheta2_2_matrix(:,i)=gradient(dsigmadtheta2_1_matrix(:,i),gtheta2(:,i));
    end
    %
    for i=1:ng1
        d2sigmadtheta12_matrix(:,i)=gradient(dsigmadtheta1_1_matrix(:,i),gtheta2(:,i));
    end
   
   dsigma=zeros(ng2,ng1);
   % now loop for each theta_i
   for j=2:ng2-1
       for i=2:ng1-1
           sigmathetatheta=[d2sigmadtheta1_2_matrix(j,i) d2sigmadtheta12_matrix(j,i);
                            d2sigmadtheta12_matrix(j,i) d2sigmadtheta2_2_matrix(j,i)];
           help = [chi_gtheta1(j,i) chi_gtheta2(j,i)];
           Term_I=sigma(j,i)*A_theta(j,i);
           Term_II=[dsigmadtheta1_1_matrix(j,i) dsigmadtheta2_1_matrix(j,i)]*B_theta(:,j,i);
           Term_III= help * sigmathetatheta*help'; Term_III=D*Term_III;
           Term_IV=C_theta(j,i);
           dsigma(j,i) = Term_I + Term_II + Term_III + Term_IV;
       end
   end
   
   dy=[];
for i=1:ng1
    dy=[dy;dsigma(:,i)];
end
   
   
   
end