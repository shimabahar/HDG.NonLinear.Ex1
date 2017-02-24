%This function calculate matrices stild
function stild = getstild()
global n ne U1 xL h Landa1 tau0

% Building matrices stild_1
stild_1 = zeros(ne-1, 1);

% Building matrices stild_2
stild_2 = zeros(ne-1, 1);

% Building matrices stild_3
stild_3 = zeros(ne-1, 1);

k = 0;
z = 0;
for i = 1:ne-1
    b = xL+(i)*h;
    g = 0;
    
    for j = 1:n 
        g = U1(j+z)*phi(j, i, b) + g;
    end
    
        stild_3(i,1) = 2*(-6*abs(Landa1(i+1)) + tau0)*g;
    k = k + n;
    z = z + n;
end;

% Building matrices stild_4
stild_4 = zeros(ne-1, 1);

for i = 1:ne-1
    stild_4(i,1) = 2*(-6*abs(Landa1(i+1)) + tau0)*Landa1(i+1);
end;

% Building matrices stild_5
stild_5 = zeros(ne-1, 1);

stild = stild_1 - stild_2 - stild_3 + stild_4 - stild_5;