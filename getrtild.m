%This function calculate matrix rtild
%rtild = -rtild_1 - rtild_2 + rtild_3
function rtild = getrtild()
global n ne U1 sigma xL h xR Say1

% Building matrix rtild_1
rtild_1 = zeros(ne, 1);

% Building matrix rtild_2
rtild_2 = zeros(ne, 1);

k = 0;
z = 0;
for i = 1:ne-1
  b = xL + i*h;
  g = 0;
    
  for j = 1:n 
    g = U1(j+z) * phi(j, i, b) + g;
  end
    
  rtild_2(i,1) = 2 * g;
  k = k + n;
  z = z + n;
end

g = 0;
for j = 1:n 
  g = U1(j + z) * phi(j, i, xR) + g;
end

rtild_2(ne, 1) = g;

% Building matrix rtild_3
rtild_3 = zeros(ne, 1);

for i = 1:ne-1
  rtild_3(i,1) = 2 * Say1(i+1);
end

rtild_3(ne, 1) = Say1(ne+1);
rtild_3 = sigma * rtild_3;

rtild = -rtild_1 - rtild_2 + rtild_3;