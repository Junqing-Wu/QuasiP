function data = QuasiP_system_identify(data, p0)
if isempty(data.l)
    data.ODE = 1;
    data.isL = false;
else
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