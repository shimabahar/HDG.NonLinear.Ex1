function out=integrateForProjectionU0(j,e,a,b)

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on per element (on [a,b])
g = 0;
for k = 1:np
  g = w(k) * U0((a+b)/2+((b-a)/2)*x(k)) * phi(j, e, (a+b)/2+((b-a)/2)*x(k)) + g;
end
out = (b-a)/2 * g;