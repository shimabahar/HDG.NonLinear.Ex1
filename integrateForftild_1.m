function out = integrateForftild_1(j,e,a,b,step)
global dt

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on per element (on [a,b])
g = 0;
for k = 1:np
  g = w(k) * f(step*dt,(a+b)/2+((b-a)/2)*x(k)) * phi(j, e, (a+b)/2+((b-a)/2)*x(k)) + g;
end
out = (b-a)/2 * g;