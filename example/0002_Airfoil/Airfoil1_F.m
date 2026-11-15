function F = Airfoil1_F(x,p,theta)
% AIRFOIL1_F  Nonlinear right-hand side of the autonomous airfoil/store model.
% Reference: Section 5.2, Eqs. (40) and (44).
% QuasiP split: A*qddot + B*qdot + C(Q)*q = F(q), with F=-K3(q).
% Only store pitch is nonlinear: F=[0;0;-10*q3^3].
% Columns of x are samples; rows 1:3 are plunge, wing pitch, store pitch.
% p and theta are accepted for callback compatibility; F does not depend on them.
% --------------------------------------------------------------------------
kbeta3 = 10;
q3 = x(3,:);
S = length(q3);
%----------------------------------------------------------
F = zeros(3,S);
F(3,:) = - kbeta3*q3.^3;
end
