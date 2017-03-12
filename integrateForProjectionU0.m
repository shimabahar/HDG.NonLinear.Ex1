function out = integrateForProjectionU0(j, e, a, b)

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on [a,b]
gl = 0;
for k = 1:np
  gl = w(k) * U0((a+b)/2 + ((b-a)/2)*x(k)) * ...
       phi(j, e, (a+b)/2 + ((b-a)/2)*x(k)) + gl;
end
out = (b-a)/2 * gl;