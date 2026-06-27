function [data, y, J] = MT_SHT_QuasiP(prob, data, u)
PHI    = data.PHI;                  % 形函数及其相关矩阵
d      = data.ode_opts.d;           % QP基频的数量
MATLABsolver = PHI.MATLABsolver;    % 1：y' = 1/Om1*g(y,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)  
                                    % 0: M*a + D*u + K*q = f(q,Om1*u,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
                                    %    M = Om1^2*A; D = Om1*B 
%----------------------------------------------------------
N      = PHI.N;                     % QP维度
NL     = PHI.nL;                    % 全体待定系数的数量
M_SHT  = PHI.M_SHT;                 % 单轨迹 积分段数
a      = data.pdim;                 % continuation paramaters [p] 的数量
Omt    = PHI.Omt;                   % 积分的基频方向，记作phi_1 = Om(Omt)*t, 周期[0,2P(Omt)pi)
%----------------------------------------------------------
switch PHI.method
    case 'HT'  % SHT_HT  超时域打靶法 
               % d==1时，只能是SHT_HT
        ys  = reshape(u(1:NL),[],N);
    case 'MF'  % SHT_MF  多频域打靶法
        Ys  = reshape(u(1:NL),[],N*M_SHT);
        ys  = reshape(PHI.Gam*Ys,[],N);
        Ys  = reshape(u(1:NL),[],N);        
end
Om  = u(1+NL:d+NL);
p   = u(d+NL+1:end);
%----------------------------------------------------------
phis  = PHI.Varphi;        % 积分段初始时刻
phi1s = PHI.Varphi1_sta;   % 积分段phi_1初始时刻
phi1e = PHI.Varphi1_end;   % 积分段phi_1终止时刻
dphi1  = PHI.dphi1;    % 定步长时积分步长 
%==========================================================
%>>>>>>>>>>>>>>>>>>>>>>>>> TI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%==========================================================
if PHI.isvar
    PSI = zeros(PHI.L_HT,N+N^2+N*d+N*a);
else
    PSI = zeros(PHI.L_HT,N);
end
if MATLABsolver
    % y' = 1/Om1*g(y,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
    parfor ii = 1 : PHI.L_HT
        ys_ii    = ys(ii,:);
        phis_ii  = phis(ii,:);
        phi1s_ii = phis_ii(Omt);
        idx_match = (phi1s == phi1s_ii);
        phi1e_ii  = phi1e(idx_match);
        [~,Psi]   = MT_SHT_MATLABsolver1(data,ys_ii,phis_ii,phi1s_ii,phi1e_ii,Om,p);
        PSI(ii,:) = Psi(end,:);
    end
else
    % M*a + D*u + K*q = f(q,Om1*u,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
    parfor ii = 1 : PHI.L_HT
        ys_ii    = ys(ii,:);
        phis_ii  = phis(ii,:);
        phi1s_ii = phis_ii(Omt);
        idx_match = (phi1s == phi1s_ii);
        phi1e_ii  = phi1e(idx_match);
        [~,Psi]   = PHI.ODEsolver(data,ys_ii(:),phis_ii(:),phi1s_ii,phi1e_ii,dphi1,Om,p);
        PSI(ii,:) = Psi(end,:);
    end
end
%==========================================================
%>>>>>>>>>>>>>>>>>>>>>>>>> BC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%==========================================================
switch PHI.method
    case 'HT'
        [y,J] = MT_SHT_HT_QuasiP(data,PSI,ys,Om);
    case 'MF'
        [y,J] = MT_SHT_MF_QuasiP(data,PSI,Ys,Om);
end
end
