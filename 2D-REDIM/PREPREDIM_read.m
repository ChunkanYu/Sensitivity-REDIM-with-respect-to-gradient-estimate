function [y0,psi_all]=PREPREDIM_read(npsi,ng1,ng2)

initial=prepredim(ng1,ng2);
ng1=size(initial,2); ng2=size(initial,1);
y0=[];
for i=1:npsi
   psi=initial(:,:,i);
   for j=1:ng1
       y0=[y0;psi(:,j)];
   end
   psi_all(:,:,i)=psi;
end



end