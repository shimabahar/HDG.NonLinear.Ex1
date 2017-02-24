% This function calculates the right hand side vector R
function R = getRhsVector(ftild,htild,ktild,rtild,stild)
global n ne 

R = zeros(3*n*ne+2*ne-1, 1);

R(1:n*ne) = ftild;
R(n*ne+1:2*n*ne) = htild;
R(2*n*ne+1:3*n*ne) = ktild;
R(3*n*ne+1:3*n*ne+ne) = rtild;
R(3*n*ne+ne+1:3*n*ne+2*ne-1) = stild;