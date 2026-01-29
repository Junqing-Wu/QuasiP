function QuasiP_arg_check(tbid, data, X_MF, OM, p0)
% QuasiP_arg_check   Basic argument checking for 'forward' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% FORWARD_ARG_CHECK(TBID, DATA, X_MF, OM, P0)
%
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% X_MF - Combination of frequency components and their amplitudes for quasi-periodic solution
% OM   - base frequencies of quasi-periodic solution
% P0   - Initial solution guess for problem parameters.

assert(numel(OM)==data.ode_opts.d, ...
  '%s: the number of OM is different from d', ...
  tbid);

assert(size(X_MF,2)==data.ode_opts.d + data.DOF, ...
  '%s: the initial value of qp is not defined', ...
  tbid);

assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);

end
