%This function calculate matrices Atild
function Atild = getAtild()
global A dt n ne Landa1 tau0 xL h U1

Atild_1 = zeros(ne*n, ne*n);
Atild_1 = (1/dt) * A;

Atild_2 = zeros(n*ne, n*ne);
Atild_2s = zeros(n, n); % sub-matrix Atild_2
k = n;
e = 1;
for l = 1:n:n*ne
  a = xL + (e-1)*h;
  b = xL + e*h;
  for i = 1:n
    for j = 1:n
      for z = 1:n
        G = integrateDphiPhi(i, j, z, e, a, b);
        Atild_2s(i, j) = Atild_2s(i, j) + G*U1(z+l-1);
      end
    end
  end
  Atild_2(l:k, l:k) = Atild_2s;
  k = k + n;
  e = e + 1;
end
Atild_2 = 6*Atild_2;

Atild_3 = zeros(n*ne, n*ne);
Atild_3s = zeros(n, n); % sub-matrix Atild_3
k = n;
z = 1;
e = 1;
for l = 1:n:n*ne 
  a = xL + (e-1)*h;
  b = xL + e*h;
  for i = 1:n
    for j = 1:n
      Atild_3s(i, j) = abs(Landa1(z+1))*phi(j, e, b)*phi(i, e, b) ...
                       + abs(Landa1(z))*phi(j, e, a)*phi(i, e, a);
    end
  end

  Atild_3(l:k, l:k) = Atild_3s;
  k = k + n;
  z = z + 1;
  e = e + 1;
end
Atild_3 = -6*Atild_3 + tau0*A;

[Atild] = Atild_1 - Atild_2 + Atild_3;