function [A,B,C,D] = cubic_spline(pS, p)
% Non-periodic cubic spline matrices
% pS: total intervals = p*S
% p: period scaling factor
% A,B,C,D: (pS+1)x(pS+1) matrices

n = pS + 1;
h = 2*p*pi / pS;

A = eye(n);

C1 = 4*h^2 * eye(n) ...
   + h^2 * diag(ones(n-1,1), 1) ...
   + h^2 * diag(ones(n-1,1), -1);

C2 = -6 * eye(n) ...
   + 3 * diag(ones(n-1,1), 1) ...
   + 3 * diag(ones(n-1,1), -1);

C1([1, n], :) = 0;
C1(1, 1) = 1;
C1(n, n) = 1;
C2([1, n], :) = 0;

C = C1 \ C2;

D1 = (-1/(3*h)) * eye(n) ...
   + (1/(3*h)) * diag(ones(n-1,1), 1);
D = D1 * C;

B1 = (-1/h) * eye(n) ...
   + (1/h) * diag(ones(n-1,1), 1);
B = B1 - h*C - h^2*D;
end