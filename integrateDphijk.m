function out=integrateDphijk(i, j,z,e,a,b)
global h

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on element 1 (on [0,h])
g = 0;
   for k = 1:np
      g = w(k) * Dphi(i, e, (h/2)*x(k)+(a+b)/2) * phi(j, e, (h/2)*x(k)+(a+b)/2)* phi(z, e, (h/2)*x(k)+(a+b)/2);
   end
out = h/2 * g;