function [sol, data] = QuasiP_read_solution(oid, run, lab, varargin)
% QUASIP_READ_SOLUTION   Read a saved QuasiP solution and toolbox data.
%
% [SOL DATA] = QUASIP_READ_SOLUTION(OID, RUN, LAB, [BRANCH_TYPE])
%
% Read the saved chart and QuasiP data for solution LAB of run RUN. If
% BRANCH_TYPE is 'BP', also recover the restart vector used for branch switching.
%
% On input:
%
% OID          : QuasiP object identifier.
% RUN          : COCO run identifier.
% LAB          : Integer solution label.
% BRANCH_TYPE  : Optional string 'BP' for branch-point restart data.
%
% On output:
%
% SOL  : Solution structure with sol.u and, for 'BP', sol.t0.
% DATA : Saved QuasiP toolbox data structure.
%
% The solution vector is reconstructed as sol.u=[C(:); Om(:); p(:)] using
% the saved index maps DATA.C_idx, DATA.OM_idx, and DATA.p_idx.
%
% See also: COCO_READ_SOLUTION, ODE_QUASIP2QUASIP, ODE_BP2QUASIP.
tbid = coco_get_id(oid, 'QuasiP');
sol.t0 = [];
if ~isempty(varargin) && strcmp(varargin{1},'BP')
    [data, chart, uidx] = coco_read_solution(tbid, run, lab, 'data', ...
        'chart', 'uidx');
    % At a BP, COCO stores the null vector used for branch switching in
    % the linear-solver chart data field lsol.v.
    cdata = coco_get_chart_data(chart, 'lsol');
    if isempty(cdata) || ~isfield(cdata, 'v') || isempty(cdata.v)
        coco_warn([], 1, 1, ...
          '%s: could not find restart data for branch-switching\n', ...
          mfilename);
    else
        sol.t0 = cdata.v(uidx);
    end
else
    [data, chart] = coco_read_solution(tbid, run, lab, 'data', 'chart');
end

% Restore the three blocks of the continuation unknown vector.
C0  = chart.x(data.C_idx);
OM0 = chart.x(data.OM_idx);
p0  = chart.x(data.p_idx); 
sol.u  = [C0(:); OM0(:); p0(:)];
end
