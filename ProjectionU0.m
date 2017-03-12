function out = ProjectionU0()
global xL h A ne n

u = zeros(n*ne, 1);

k = 0;
for i = 1:ne
    a = xL + (i-1)*h;
    b = xL + i*h;
    for j = 1:n
        u(k+j) = integrateForProjectionU0(j, i, a, b);
    end
    k = k+n;
end

out = A \ u;