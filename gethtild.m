%This function calculate matrix htild
%htild = -htild_1 - htild_2 + htild_3 + htild_4
function htild = gethtild()
global n ne Q1 P1 A B Say1 sigma xL h

% Building matrix htild_1
htild_1 = A * P1;

% Building matrix htild_2
htild_2 = B * Q1;

% Building matrix htild_3
htild_3 = sigma * A * Q1;

% Building matrix htild_4
htild_4 = zeros(n*ne, 1);

k=0;
z = 0;
for i = 1:ne
  a = xL + (i-1)*h;
  b = xL + i*h;
   
  for j = 1:1:n
    htild_4(k+j,1) = (1 -sigma)*Say1(i + 1)*phi(j, i, b) + (-1 - sigma)*Say1(i)*phi(j, i, a);
  end
    
  k = k + n;
  z = z + n;
end

htild = -htild_1 - htild_2 + htild_3 + htild_4;