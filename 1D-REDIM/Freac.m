function F=Freac(psi,k)
  psi1=psi(1); psi2=psi(2); psi3=psi(3);
  k1=k(1); k2=k(2); k3=k(3);
  F(1,1)=-k1*psi1;
  F(1,2)=k1*psi1-k2*psi2;
  F(1,3)=k2*psi2-k3*psi3;
end