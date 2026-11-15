function sol  = QuasiP_init_sol(data,Xk_MF,OM,p0)
%QUASIP_INIT_SOL  Convert Starter data to the COCO unknown vector.
%
% INITIAL_C_byPHI maps the multi-frequency Fourier representation Xk_MF to
% the coefficient ordering C used by the selected PHI method. The solution
% vector then follows the NAE ordering u=[C(:); Om(:); p0(:)].
%
% On input:
%   data  : Initialized QuasiP data and PHI discretization.
%   Xk_MF : Harmonic combinations k_tilde and Fourier coefficients X_c(k_tilde).
%   OM    : Initial omega_1^d.
%   p0    : Initial model parameters.
%
% On output:
%   sol   : Structure with sol.u ready for coco_add_func.
C = INITIAL_C_byPHI(data,Xk_MF,OM);

sol.u  = [C(:); OM(:); p0(:)];
end