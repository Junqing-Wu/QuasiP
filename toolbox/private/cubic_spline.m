function [A,B,C,D] = cubic_spline(pS, p)
%CUBIC_SPLINE  Construct cubic-spline coefficient matrices on a periodic mesh.
%
% Build the interval-wise cubic interpolation coefficients used by the finite
% difference and post-processing routines. The period multiplier p determines
% the repeated interval layout and the cyclic boundary equations.
%
% On input:
%   pS : Number of mesh intervals in the periodic grid.
%   p  : Period multiplier; the interval spans 2*p*pi.
%
% On output:
%   A : Cubic-spline value coefficients for each interval.
%   B : Cubic-spline first-order coefficients.
%   C : Cubic-spline second-order coefficients.
%   D : Cubic-spline third-order coefficients.
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