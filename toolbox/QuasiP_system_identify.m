function data = QuasiP_system_identify(data, p0)
%QUASIP_SYSTEM_IDENTIFY  Identify the ODE representation used by Generator.
%
% If l is empty, QuasiP uses a first-order model. Otherwise l(p0) reports
% which A(p), B(p), and C(p) terms are present in
% l=A(p)*ddx+B(p)*dx+C(p)*x and which depend on p. The flags created here
% select the correct residual, Jacobian, amplitude, and stability paths.
%
% On input:
%   data : Model callbacks f, fx, fp, l and ode_opts.
%   p0   : Parameter vector used to evaluate l(p0).
%
% On output:
%   data : Updated ODE-order and callback-availability flags.
if isempty(data.l)
    data.ODE = 1;
    data.isL = false;
else
    % Query which second-order matrices and parameter derivatives are active.
    [isABC,isABCwithp] = data.l(p0);
    if isABC(1) == 1
        data.ODE = 2;
    else
        data.ODE = 1;
    end
    data.isL        = true;
    data.isABC      = isABC;
    data.isABCwithp = isABCwithp;
end
%----------------------------------------------------------
data.isp      = ~isempty(p0);
data.isfwithp = ~isempty(data.fp);
data.isfwithu = data.ode_opts.isfwithu;
end