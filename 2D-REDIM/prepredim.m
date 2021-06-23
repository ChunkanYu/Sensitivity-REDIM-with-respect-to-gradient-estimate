function [initial]=prepredim(ng1,ng2)

npsi=3;

upper_original=importdata('./prepredim/inipro_chi_0p6.mat');upper_original=reshape(upper_original,[],npsi)';

below_original=importdata('./prepredim/inipro_chi_1000.mat');below_original=reshape(below_original,[],npsi)';


ptheta1=1; ptheta2=3;

theta1_original=upper_original(ptheta1,:); 
theta1_min=min(theta1_original)*1.005; theta1_max=0.995*max(theta1_original);
theta1=linspace(theta1_min,theta1_max,ng1);

for i=1:npsi
    upper(i,:) = interp1(theta1_original,upper_original(i,:),theta1);
    below(i,:) = interp1(theta1_original,below_original(i,:),theta1);
end

for j=1:npsi
for i=1:ng1
    initial(:,i,j)=linspace(below(j,i),upper(j,i),ng2);
end
end

end
