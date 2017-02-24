%This function calculate matrices C1tild
function [C1tild] = getC1tild()
global n ne Landa1 X tau0 U

C1tild_1 = zeros(ne*n, ne+1);
C1tild_1s = zeros(n, 2); % sub-matrix C1tild_1
e = 1;
for s = 1:n:n*ne
  for i = 1:n
     for j = 1:2
     C1tild_1s(i, j) = (-1)^j* Landa1(e+j-1,1)*phi(i, e, X(e+j-1));
     end
  end
  C1tild_1(s:s+n-1, e:e+1) = C1tild_1s;
  e = e+1;
end
C1tild_1 = 6*C1tild_1(:, 2:ne);

C1tild_21 = zeros(ne*n, ne+1);
C1tild_21s = zeros(n, 2); % sub-matrix C1tild_21
e = 1;
kk = 1;
g = 0;
for s = 1:n:n*ne
  V = U(kk:kk+n-1,1);
  kk = kk+n;
  for i = 1:n
      jj = 1;
     for j = 1:2
         if Landa1(j+e-1,1) >= 0
             jj = 0;
         end
         for k = 1:n
             g = g + phi(k,e,X(e+j-1))*V(k,1);
         end;
         % b khatere moshtagh 6 zarb shode
     C1tild_21s(i, j) = (-1)^jj*6*g*phi(i, e, X(e+j-1));
     g = 0;
     end
  end
  C1tild_21(s:s+n-1, e:e+1) = C1tild_21s;
  e = e+1;
end
C1tild_21 = C1tild_21(:, 2:ne);

C1tild_22 = zeros(ne*n, ne+1);
C1tild_22s = zeros(n, 2); % sub-matrix C1tild_22
e = 1;
for s = 1:n:n*ne
  for i = 1:n
      jj = 1;
     for j = 1:2
         if Landa1(j+e-1,1)>=0
             jj=0;
         end
     C1tild_22s(i, j) =Landa1(j+e-1,1)*(-1)^jj*6*phi(i, e, X(e+j-1));
     end
  end
  C1tild_22(s:s+n-1, e:e+1) = C1tild_22s;
  e = e+1;
end
C1tild_22 = C1tild_22(:, 2:ne);

C1tild_23 = zeros(ne*n, ne+1);
C1tild_23s = zeros(n, 2); % sub-matrix C1tild_23
e = 1;
for s = 1:n:n*ne
  for i = 1:n
     for j = 1:2
     C1tild_23s(i, j) =(-6*abs(Landa1(j+e-1,1))+tau0)*(phi(i, e, X(e+j-1)));
     end
  end
  C1tild_23(s:s+n-1, e:e+1) = C1tild_23s;
  e = e+1;
end
C1tild_23 = C1tild_23(:, 2:ne);

C1tild_2 = C1tild_21 - C1tild_22 - C1tild_23;
C1tild = C1tild_1 + C1tild_2;