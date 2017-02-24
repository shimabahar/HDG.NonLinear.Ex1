%This function calculate matrices E2tild
function [E2tild] = getE2tild()
global ne Landa1 tau0 n U X

E2tild_1=zeros(ne-1,ne-1);

E2tild_21 = zeros(ne-1, ne-1);
e = 1;
k=1;
kk=n;
g=0;
gg=0;
for s = 1:ne-1
  V=U(k:k+n-1,1);
  VV=U(kk+1:kk+n,1);
  k=k+n;
  kk=kk+n;
  if Landa1(s+1)>=0
     jj=0;
  end
  for i=1:n
     g=g+phi(i,e,X(s+1))*V(i,1);
  end;
  for i=1:n
     gg=gg+phi(i,e+1,X(s+1))*VV(i,1);
  end;
  E2tild_21s =(-1)^jj*6*(g+gg);
  E2tild_21(s,s) = E2tild_21s;
  e = e+1;
  g=0;
  gg=0;
end

E2tild_22 = zeros(ne-1, ne-1);
for s = 1:ne-1
   jj=1;
   if Landa1(s+1)>=0
      jj=0;
   end
  E2tild_22s =Landa1(s+1)*(-1)^jj*6*2;
  E2tild_22(s,s) = E2tild_22s;
end

E2tild_23 = zeros(ne-1, ne-1);
for s = 1:ne-1
     C1tild_23s =(-6*abs(Landa1(s+1))+tau0)*2;
  E2tild_23(s,s) = C1tild_23s;
end

E2tild_2=E2tild_21-E2tild_22-E2tild_23;
E2tild=E2tild_1+E2tild_2;