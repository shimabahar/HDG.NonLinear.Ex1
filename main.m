%{
Matlab code to solve KdV equation using HDG method.

More details can be found in the following paper:
A Hybridized Discontinuous Galerkin Method for the Nonlinear
Korteweg–de Vries Equation - Ali Samii et al.

Author(s):
Date:
%}

clc
clear all

tic

global n ne xL xR dt T tau0 sigma h X Landa U U1 Q1 P1 step A B Landa1 Say1

%%%%%% Problem configuration %%%%%%
degree = 1; % max degree of the base functions (phi)
n = degree + 1;
ne = 10; % number of elements in x domain
xL = 0;  % left hand side
xR = pi; % right hand side
dt = 10^-4;
T = 0.1;
tau0 = 20;
sigma = 1;
h = double((xR-xL)/ne);
X = xL:h:xR;
tol = 10^-3; % error tolerance

%%%% Build coefficient matrices and the left hand side matrix %%%%
[A B] = getAandB();
[D2] = getD2();
[C2 C3] = getC2andC3();
[C4T C5T C6T] = getC4TandC5TandC6T();
[E1] = getE1();

U = ProjectionU0();
Q = ProjectionQ0();
P = ProjectionP0();

Landa = ones(ne+1, 1);
Landa(1) = guxL(0);
Landa(ne+1) = guxR(0);
Say = ones(ne+1,1);
Say(1) = gqxL(0);

%M = getLhsMatrix(A, B, D1, D2, C1, C2, C3, C4T, C5T, C6T, E1, E2);

U1 = U;
Q1 = Q;
P1 = P;
Landa1 = Landa;
Say1 = Say;

%%%% Solve the system and calculate the error  %%%%
for step = 0:T/dt 
  for j = 1:1000
    Atild = getAtild();
    C1tild = getC1tild();
    E2tild = getE2tild();

    ftild = getftild();
    htild = gethtild();
    ktild = getktild();
    rtild = getrtild();
    stild = getstild();

    R = getRhsVector(ftild,htild,ktild,rtild,stild);
    M = getLhsMatrix(Atild, C1tild, A, B, D2, C2, C3, C4T, C5T, C6T, E1, E2tild);
    Answer = M \ R;

    U1 = U1 + Answer(1:n*ne);
    Q1 = Q1 + Answer(n*ne+1:2*n*ne);
    P1 = P1 + Answer(2*n*ne+1:3*n*ne);
    landa1(2:ne) = Landa1(2:ne) + Answer(3*n*ne+1:3*n*ne+ne-1);
    Say1(2:ne+1) = Say1(2:ne+1) + Answer(3*n*ne+ne:3*n*ne+2*ne-1);
    
    if (U1 <= tol) & (Q1 <= tol) & (P1 <= tol) & (Landa1 <= tol) & (Say1 <= tol)
      break;
    end
  end

  U = U1;
  Q = Q1;
  P = P1;
  Landa = Landa1;
  Say = say1;

  Landa(1) = guxL(step*dt);
  Landa(ne+1) = guxR(step*dt);
  Say(1) = gqxL(step*dt);

  error = getL2Error(U0, step);
  disp(error);
end

toc