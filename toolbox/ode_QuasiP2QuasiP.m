function prob = ode_QuasiP2QuasiP(prob,oid,varargin)
% ODE_QUASIP2QUASIP   Restart QuasiP continuation from a saved solution.
%
% PROB = ODE_QUASIP2QUASIP(PROB, OID, RUN, [SOID], LAB)
%
% Read a previously saved QuasiP solution and use it as the initial point
% of a new continuation problem.
%
% On input:
%
% PROB : COCO continuation problem structure.
% OID  : Target QuasiP object identifier.
% RUN  : Source continuation run.
% SOID : Optional source object identifier.
% LAB  : Saved solution label.
%
% On output:
%
% PROB : COCO problem with the saved QuasiP solution restored.
tbid = coco_get_id(oid, 'QuasiP'); % Create toolbox instance identifier
str  = coco_stream(varargin{:});   % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
    soid = str.get;
else
    soid = oid;
end
lab = str.get;
[sol,data] = QuasiP_read_solution(soid, run, lab);
data = QuasiP_pnames_update(tbid, data, data.pnames_np);
QuasiP_print_equation_table(data);
prob = ode_construct_QuasiP(prob, tbid, data, sol);
end