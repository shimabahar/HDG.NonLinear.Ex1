%This function calculate matrix ktild
%ktild = -ktild_1 - ktild_2 + ktild_3
function ktild = getktild()
global n ne U1 Q1 A B Landa1 xL h

% Building matrix ktild_1
ktild_1 = A * Q1;

% Building matrix ktild_2
ktild_2 = B * U1;

% Building matrix ktild_3
ktild_3 = zeros(n*ne, 1);

k = 0;
z = 0;
for i = 1:ne
  a = xL + (i-1)*h;
  b = xL + i*h;
    
  for j = 1:n
    htild_4(k+j, 1) = Landa1(i+1) * phi(j, i, b) - Landa1(i) * phi(j, i, a);
  end
    
  k = k + n;
  z = z + n;
end

ktild = -ktild_1 - ktild_2 + ktild_3;