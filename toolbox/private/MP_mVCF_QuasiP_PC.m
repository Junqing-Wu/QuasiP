function [data, y, J] = MP_mVCF_QuasiP_PC(prob, data, u)
PHI    = data.PHI;          % 形函数及其相关矩阵
d      = data.ode_opts.d;   % QP基频的数量
e      = data.ode_opts.e;   % 激励基频的数量
nL     = PHI.nL;            % 与QP相关的待定系数的数量
a      = data.pdim;         % continuation paramaters [p] 的数量
%--------------------------------------------------------------------------
y1   = zeros(d-e,1);        % 未知基频的相位条件PC的残差
Jx1  = zeros(d-e,nL);       % PC关于Xd的导数  
Jom1 = zeros(d-e,d);        % PC关于Om的导数 
Jp1  = zeros(d-e,a);        % PC关于p 的导数 
%==========================================================
Xd   = u(1:nL);             % nL*1
if d>e
    kk = e+1;
    for ii = 1:d-e
        if data.ODE == 1
            Jx1(ii,:) = transpose(Xd)*transpose(PHI.Lam1{kk});
        else
            Jx1(ii,:) = transpose(Xd)*...
                        ( transpose(PHI.Lam1{kk}) + ...
                          transpose(PHI.Lam2{kk})*PHI.Lam1{kk} );
        end
        kk = kk+1;
    end
end
%==========================================================
y = y1;
J = [Jx1, Jom1, Jp1];
end