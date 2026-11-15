function DFDX = Airfoil1_DFDX(x,p,theta)
% AIRFOIL1_DFDX  State Jacobian of the nonlinear right-hand side.
% Reference: Section 5.2, Eq. (44), with F=-K3(q).
% The only nonzero derivative is dF3/dq3=-30*q3^2.
% Return a 3-by-3-by-S array for S input samples.
% This callback differentiates F only; the linear matrices are handled separately.
% --------------------------------------------------------------------------
kbeta3 = 10;
q3 = x(3,:);
S = length(q3);
%----------------------------------------------------------
DFDX = zeros(3,3,S);
DFDX(3,3,:) = - 3 * kbeta3 * q3.^2;
%----------------------------------------------------------
end