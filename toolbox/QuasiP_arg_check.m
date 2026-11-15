function QuasiP_arg_check(tbid, data, Xk_MF, OM, p0)
%QUASIP_ARG_CHECK  Validate the dimensions and consistency of initial data.
%
% This function is called before COCO registration. It prevents a mismatch
% between the requested number of base frequencies, the state dimension, the
% harmonic-index matrix, the coefficient matrix, and the model-parameter names.
%
% On input:
%   tbid : QuasiP identifier used in error messages.
%   data : Settings and model metadata.
%   Xk_MF: Initial multi-frequency representation with k_wave and Xc.
%   OM   : Initial vector of d base frequencies.
%   p0   : Initial model-parameter vector.
%
% On output:
%   None. An assertion error stops construction when a check fails.
assert(numel(OM)==data.ode_opts.d, ...
  '%s: the number of OM is different from d', ...
  tbid);

% Each harmonic-index row must have one entry per base frequency.
assert(size(Xk_MF.k_wave,2)==data.ode_opts.d, ...
  '%s: the number of harmonics in k_wave is different from d', ...
  tbid);

% The initial coefficient matrix must have one column per state component.
assert(size(Xk_MF.Xc,2)==data.ode_opts.n, ...
  '%s: the number of amplitude columns in Xc is different from n', ...
  tbid);

assert(size(Xk_MF.Xc,1)==size(Xk_MF.k_wave,1), ...
  '%s: Xc and k_wave have inconsistent numbers of rows', ...
  tbid);

assert(numel(p0)==numel(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);

end