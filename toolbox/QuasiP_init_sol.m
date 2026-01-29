function sol  = QuasiP_init_sol(data,X_MF,OM,p0)
% QUASIP_INIT_SOL   Build initial solution guess.
%
% Construct initial solution guess for 'QuasiP' toolbox.
%
% SOL = QUASIP_INIT_SOL(DATA, X_MF, OM, P0)
%
% DATA - Toolbox data structure.
% X_MF - Combination of frequency components and their amplitudes for quasi-periodic solution
% OM   - base frequencies of quasi-periodic solution
% P0   - Initial solution guess for problem parameters.

C = INITIAL_C_byPHI(data,X_MF,OM);

sol.u  = [C(:); OM(:); p0(:)];
end