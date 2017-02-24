% This function calculate matrices C1, C2, and C3
function [C2 C3] = getC2andC3()
global n ne sigma X

CL2 = zeros(ne*n, ne+1);
CL3 = zeros(ne*n, ne+1);
CLs = zeros(n, 2); % sub-matrix CL (C Large)

for i = 1:n
  for j = 1:2
    CLs(i, j) = phi(i, 2, X(j+1));
  end
end

CLs2(:, 1) = (-1-sigma) * CLs(:, 1);
CLs2(:, 2) = (1-sigma) * CLs(:, 2);
CLs3(:, 1) = -1 * CLs(:, 1);
CLs3(:, 2) = +1 * CLs(:, 2);

p = 1;
for s = 1:n:n*ne
  CL2(s:s+n-1, p:p+1) = CLs2;
  CL3(s:s+n-1, p:p+1) = CLs3;
  p = p+1;
end

% Building matrices Ci's by eliminating extra columns from CLi's
C2 = CL2(:, 2:ne+1);
C3 = CL3(:, 2:ne);