function prob = ode_BP2QuasiP(prob,oid,varargin)
% ODE_BP2QUASIP   Restart QuasiP continuation from a saved branch point.
%
% PROB = ODE_BP2QUASIP(PROB, OID, RUN, [SOID], LAB)
%
% Restore the saved QuasiP solution and the secondary tangent stored in lsol.v.
%
% On input:
%
% PROB : COCO continuation problem structure.
% OID  : Target QuasiP object identifier.
% RUN  : Source continuation run.
% SOID : Optional source object identifier.
% LAB  : Saved branch-point label.
%
% On output:
%
% PROB : COCO problem with the branch-point solution restored.
%
% See also: QUASIP_READ_SOLUTION, ODE_QUASIP2QUASIP.
tbid = coco_get_id(oid, 'QuasiP'); % Create toolbox instance identifier
str  = coco_stream(varargin{:});   % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
    soid = str.get;
else
    soid = oid;
end
lab = str.get;
[sol,data] = QuasiP_read_solution(soid, run, lab, 'BP');
if isempty(sol.t0) || any(~isfinite(sol.t0)) || norm(sol.t0) == 0
    error('ode_BP2QuasiP:InvalidRestartTangent', ...
        'No usable BP restart tangent was saved for label %d.', lab);
end
% lsol.v already points along the secondary branch. Refine and follow it.
prob = coco_set(prob, 'cont', 'branch', 'no switch');

if ~strcmp(oid, soid)
   pnames = data.pnames;
   data = QuasiP_pnames_update(tbid, data, pnames);
end
QuasiP_print_equation_table(data);
prob = ode_construct_QuasiP(prob, tbid, data, sol);
end
