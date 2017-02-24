%This function calculate matrices ftild
%ftild = ftild_1 - ftild_2 + ftild_3 + ftild_4 - ftild_5
function ftild = getftild()
global n ne step xL xR h B U U1 dt A Landa1 tau0

% Building matrices ftild_1
ftild_1 = zeros(n*ne,1);

k = 0;
for i = 1:ne
    a = xL+(i-1)*h;
    b = xL+(i)*h;
    
    for j = 1:1:n
        ftild_1(k+j,1) = integrateForftild_1(j,i,a,b,step);
    end;
    k = k+n;
end;

% Building matrices ftild_2
ftild_2 = zeros(n*ne,1);

ftild_2=(1/dt)*A*(U1-U);

% Building matrices ftild_3
ftild_3 = zeros(n*ne,1);

ftild_3=B'*(U1);

% Building matrices ftild_4
ftild_4 = zeros(n*ne,1);

k=0;
z = 0;
for i=1:ne
    a = xL+(i-1)*h;
    b = xL+(i)*h;
    
    u = zeros(n, 1);
    for j = 1:n  % differenet base functions
        u(j) = U1(j+z);
    end
    
    for j = 1:1:n
        ftild_1(k+j,1) = integrateForftild_4(j,i,a,b,u);
    end;
    k = k+n;
    z=z+n;
end;

% Building matrices ftild_5
%ftild_5 = ftild_51 + ftild_52 - ftild_53
ftild_51 = zeros(n*ne,1);
ftild_52 = zeros(n*ne,1);
ftild_53 = zeros(n*ne,1);

%calculate matrices ftild_51 & ftild_53
k=0;
z = 0;
for i=1:ne
    a = xL + (i-1)*h;
    b = xL + (i)*h;
    
    for j = 1:1:n
        ftild_51(k+j,1) = Landa1(i+1)*phi(j, i, b) + Landa1(i)*phi(j, i, a);
        ftild_53(k+j,1) = (-6*abs(Landa1(i+1)) + tau0)*(Landa1(i+1)*phi(j, i, b)) - (Landa1(i)*phi(j, i, a))*(-6*abs(Landa1(i)) + tau0);
    end;
    
    k = k + n;
    z = z + n;
end;

%calculate matrices ftild_51 & ftild_53
k=0;
z = 0;
for i=1:ne
    a = xL+(i-1)*h;
    b = xL+(i)*h;
    g1 = 0;
    g2 = 0;
    
    for j = 1:n 
        g1 = U1(j+z)*phi(j,i,a) + g1;
        g2 = U1(j+z)*phi(j,i,b) + g2;
    end
    
    for l = 1:n
        ftild_52(k+l,1) = (-6*abs(Landa1(i+1)) + tau0)*g2*phi(j, i, b) - (-6*abs(Landa1(i)) + tau0)*g1*phi(j, i, a);
    end;
    k = k + n;
    z = z + n;
end;

ftild_5 = ftild_51 + ftild_52 - ftild_53;
ftild = ftild_1 - ftild_2 + ftild_3 + ftild_4 - ftild_5;

