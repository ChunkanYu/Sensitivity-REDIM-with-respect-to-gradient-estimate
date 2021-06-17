function dpsidtheta_perp=gram_schmidt(dpsidtheta_vector)

ns=size(dpsidtheta_vector,1);

xx=eye(ns); xxx=xx(:,1:end-1); 

vv=[dpsidtheta_vector,xxx];

u=zeros(ns,ns);

u(:,1:1)=vv(:,1:1);

for k=2:ns
    v=vv(:,k);
    yyy=v*0;
    for i=1:k-1
       yyy=yyy+proj(u(:,i),v); 
    end
    u(:,k)=v-yyy;
    dpsidtheta_perp(:,k-1)=u(:,k)/norm(u(:,k));
end


end


function proj_uv=proj(u,v)
proj_uv = ((u'*v)/(u'*u))*u;
end