function Cons = MT_SHT_get_MDK_Cons(data, p, n, a, Om1)
%MT_SHT_GET_MDK_CONS  Precompute constant data for the SHT integrator.
%
% Evaluate the linear callback l at p and organize A, B, C and their
% parameter derivatives in a compact structure. MT_SHT_MDKF_var reuses this
% data at every integration point.
%
% On input:
%   data : QuasiP model callbacks and ODE flags.
%   p    : Model-parameter vector.
%   n    : Original ODE state dimension.
%   a    : Number of model parameters.
%   Om1  : Selected integration frequency.
%
% On output:
%   Cons : Precomputed mass, damping, stiffness, and derivative data used by
%          MT_SHT_MDKF_var.

[~,~,A,B,C,Ap,Bp,Cp] = data.l(p);
%-------------------------------------------------------------------------------------
if data.isABC(2) ~= 1
    B = zeros(n,n);end
if data.isABC(3) ~= 1
    C = zeros(n,n);end
if data.isABCwithp(1) ~= 1
    Ap = zeros(n,n,a);end
if data.isABCwithp(2) ~= 1
    Bp = zeros(n,n,a);end
if data.isABCwithp(3) ~= 1
    Cp = zeros(n,n,a);end
%-------------------------------------------------------------------------------------
Cons.A         = A;
Cons.Aw        = A*Om1^2;
Cons.invA      = inv(A);
%-------------------------------------------------------------------------------------
Cons.B         = B;
Cons.Bw        = B*Om1;
%-------------------------------------------------------------------------------------
Cons.C         = C;
%-------------------------------------------------------------------------------------
if data.isp
    Cons.A_p     = Ap;
    Cons.B_p     = Bp;
    Cons.C_p     = Cp;
else
    Cons.A_p     = zeros(n,n,a);
    Cons.B_p     = zeros(n,n,a);
    Cons.C_p     = zeros(n,n,a);
end
end