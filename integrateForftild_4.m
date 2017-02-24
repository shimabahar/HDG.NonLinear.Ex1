function out=integrateForftild_4(j,e,a,b,u)
global n dt

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on element 1 (on [0,h])
gl = 0;

for i = 1:np
  Utild = 0;
  for k = 1:n
    Utild = u(k) * phi(k, e, (a+b)/2+(b-a)/2*x(i)) + Utild;
  end
  gl = w(i) * (Utild^2 *phi(j, e, (a+b)/2+(b-a)/2*x(i))) + gl;
end

out = ((b-a)/2)^2 * gl;