function [data, y, J] = MP_mVCF_QuasiP_PC(prob, data, u)
%MP_MVCF_QUASIP_PC  Evaluate orthogonal phase conditions for m-VCF.
%
% When d>e, the response frequencies Om(e+1:d) are unknown. Each unknown
% frequency introduces one phase condition because a torus can be shifted in
% the corresponding hyper-time direction without changing the physical torus.
% The code forms one row for each such frequency using PHI.Lam1 and, for
% second-order systems, PHI.Lam2.
%
% On input:
%   prob : COCO problem structure.
%   data : Model and PHI data.
%   u    : Unknown vector whose first nL entries are the m-VCF coefficients.
%
% On output:
%   data : Unchanged or updated callback data.
%   y    : Orthogonal phase-condition residual.
%   J    : Derivative of the phase conditions with respect to u.
PHI    = data.PHI;          % PHI basis, transform, and differentiation data
d      = data.ode_opts.d;   % Number d of response base frequencies
e      = data.ode_opts.e;   % Number e of excitation frequencies
nL     = PHI.nL;            % Total coefficient count nL=N*L
a      = data.pdim;         % Number a of model parameters p
pos    = PHI.pos;
%--------------------------------------------------------------------------
R1   = zeros(d-e,1);        % OPC residual for the d-e unknown frequencies
JXd1 = zeros(d-e,nL);       % dOPC/dC  
Jom1 = zeros(d-e,d);        % dOPC/dOm 
Jp1  = zeros(d-e,a);        % dOPC/dp 
%==========================================================
Xd   = u(1:nL);             % nL*1
if d>e
    kk = e+1;
    for ii = 1:d-e
        if data.ODE == 1
            JXd1(ii,:) = 1e5* transpose(Xd)*transpose(PHI.Lam1{pos(kk)});
        else
            JXd1(ii,:) = 1e5* transpose(Xd)*...
                        ( transpose(PHI.Lam1{pos(kk)}) + ...
                          transpose(PHI.Lam2{pos(kk)})*PHI.Lam1{pos(kk)} );
        end
        kk = kk+1;
    end
end
%==========================================================
y = R1;
J = [JXd1, Jom1, Jp1];
end