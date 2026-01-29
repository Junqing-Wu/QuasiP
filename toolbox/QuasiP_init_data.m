function data = QuasiP_init_data(data, OM, p0)
% QuasiP_init_data   Initialize toolbox data for an instance of 'QuasiP'.
%
% Populate remaining fields of the toolbox data structure used by 'QuasiP'
% function objects.
%
% DATA = QUASIP_INIT_DATA(DATA, OM, P0)
%
% DATA - Toolbox data structure.
% OM   - base frequencies of quasi-periodic solution.
% P0   - Initial solution guess for problem parameters.

data.cdim   = data.nL;   % Number of undetermined coefficients of qp 
data.omdim  = numel(OM); % Number of base frequencies 
data.pdim   = numel(p0);  % Number of problem parameters

data.C_idx  = 1:data.cdim;
data.OM_idx  = data.C_idx(end) + (1:data.omdim);
data.p_idx  = data.OM_idx(end) + (1:data.pdim);

% data = coco_func_data(data);

end