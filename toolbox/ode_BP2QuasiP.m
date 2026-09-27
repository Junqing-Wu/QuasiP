function prob = ode_BP2QuasiP(prob,oid,varargin)
% ODE_BP2QuasiP   Append 'QuasiP' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = ODE_BP2QUASIP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB  - Continuation problem structure.
% OID   - Target object instance identifier (string).
% RUN   - Run identifier (string).
% SOID  - Source object instance identifier (string).
% LAB   - Solution label (integer).

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

if ~strcmp(oid, soid)
   pnames = data.pnames;
   data = QuasiP_pnames_update(tbid, data, pnames);
end
QuasiP_print_equation_table(data);
prob = ode_construct_QuasiP(prob, tbid, data, sol);
end