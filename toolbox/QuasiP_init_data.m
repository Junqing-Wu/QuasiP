function data = QuasiP_init_data(data, OM, p0)
%QUASIP_INIT_DATA  Define dimensions and index maps for the NAE variables.
%
% PHI.nL is L*N, where N is the number of undetermined system variables and
% L is the number of retained coefficients per variable. The unknown vector is
% u=[C(:); Om(:); p(:)]. This routine records the block dimensions and the
% index maps used by residuals, OPCs, parameter monitors, and saved charts.
%
% On input:
%   data : Initialized QuasiP data containing PHI.nL, d, and e.
%   OM   : Initial omega_1^d vector.
%   p0   : Initial model-parameter vector.
%
% On output:
%   data : Updated cdim, omdim, pdim, C_idx, OM_idx, p_idx, kOM_idx, and
%          uOM_idx fields.
d = data.ode_opts.d;
e = data.ode_opts.e;

data.cdim   = data.PHI.nL;   % Number of undetermined coefficients of qp 
data.omdim  = numel(OM);     % Number of base frequencies 
data.pdim   = numel(p0);     % Number of problem parameters

% Define contiguous blocks in u = [C; Om; p].
data.C_idx   = 1:data.cdim;
data.OM_idx  = data.C_idx(end) + (1:data.omdim);
data.p_idx   = data.OM_idx(end) + (1:data.pdim);
% Split frequencies into prescribed excitation frequencies and unknown response frequencies.
data.kOM_idx = data.OM_idx(1:e);
data.uOM_idx = data.OM_idx(1+e:d);
end