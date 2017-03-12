function out = integrateDphiPhi(i, j, z, e, a, b)
global h

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on [a, b]
gl = 0;
for k = 1:np
  gl = w(k) * Dphi(i, e, (h/2)*x(k)+(a+b)/2) ...
       * phi(j, e, (h/2)*x(k)+(a+b)/2)* phi(z, e, (h/2)*x(k)+(a+b)/2) ...
       + gl;
end
out = h/2 * gl;