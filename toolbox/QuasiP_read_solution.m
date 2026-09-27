function [sol, data] = QuasiP_read_solution(oid, run, lab, varargin)
% QUASIP_READ_SOLUTION   Read 'QuasiP' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'QuasiP' toolbox instance
% identifier from solution file and construct solution structure.
%
% [SOL DATA] = QUASIP_READ_SOLUTION(OID, RUN, LAB, VARARGIN)
% VARARGIN = { [BRANCH_TYPE] }
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).
% BRANCH_TYPE - Branch type ('BP').

tbid = coco_get_id(oid, 'QuasiP');
sol.t0 = [];
if ~isempty(varargin) && strcmp(varargin{1},'BP')
    [data, chart, uidx] = coco_read_solution(tbid, run, lab, 'data', ...
        'chart', 'uidx');
    cdata = coco_get_chart_data(chart, 'lsol');
    if isempty(cdata)
      coco_warn([], 1, 1, ...
        'Could not find restart data for branch-switching\n');
    else
      sol.t0 = cdata.v(uidx);
    end
else
    [data, chart] = coco_read_solution(tbid, run, lab, 'data', 'chart');
end

C0  = chart.x(data.C_idx);
OM0 = chart.x(data.OM_idx);
p0  = chart.x(data.p_idx); 
sol.u  = [C0(:); OM0(:); p0(:)];
end
