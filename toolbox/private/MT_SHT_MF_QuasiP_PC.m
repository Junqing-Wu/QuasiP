function [data, y, J] = MT_SHT_MF_QuasiP_PC(prob, data, u)
% MF-OPC relation:
%   OPC_i = <phase-direction derivative, continuation correction> = 0.
% The code evaluates this row using the initial multi-frequency coefficients.
%MT_SHT_MF_QUASIP_PC  Build phase conditions for multi-frequency shooting.
%
% The unknown response frequencies require one phase condition per free
% hyper-time direction. This routine evaluates the selected phase-direction
% derivative on the initial multi-frequency coefficient blocks Ys and returns
% the rows that are appended to the shooting continuation problem.
%
% On input:
%   prob : COCO problem structure.
%   data : QuasiP model and PHI data.
%   u    : Shooting unknown vector.
%
% On output:
%   data : Updated callback data.
%   y    : Multi-frequency shooting phase-condition residual.
%   J    : Its Jacobian with respect to the active unknowns.
PHI   = data.PHI;
NL    = PHI.nL;
N     = PHI.N;
d     = data.ode_opts.d;
e     = data.ode_opts.e;
a     = data.pdim;
%----------------------------------------------------------
R1   = zeros(d-e,1);        % OPC residual for unknown response frequencies
JY1  = zeros(d-e,NL);       % dOPC/dC
Jom1 = zeros(d-e,d);        % dOPC/dOm
Jp1  = zeros(d-e,a);        % dOPC/dp
%----------------------------------------------------------
M_SHT = PHI.M;
Ys  = reshape(u(1:NL),[],N*M_SHT);
ys  = reshape(PHI.Gam*Ys,[],N);
Ys  = reshape(u(1:NL),[],N);
Om  = u(1+NL:d+NL);
p   = u(d+NL+1:end);
Om1 = Om(PHI.Omt);
%--------------------------------------------------------------------------------------------------------------------
Omt     = PHI.Omt;
Omnot   = PHI.Omnot;
Rho     = Om(Omnot)./Om(Omt);
Varphi  = PHI.Varphi;
if d == 1
    Theta   = Varphi;
else
    Rhophi1 = Varphi(:,Omt).*Rho';
    Theta   = Varphi;
    Theta(:,Omnot) = Theta(:,Omnot) + Rhophi1;
end
%============================================================
%>>>>>>>>>>>>>>>                dYs1 orthogonal to Ys1r <<<<<<<<<<<<<<
%>>>>>>>>>>>>>>> (I_N \otimes Xi)Ys1 orthogonal to Ys1r <<<<<<<<<<<<<<
%============================================================
if d > e
    for mm = 1:PHI.M
        kk = e + 1;
        loca0   = PHI.local_HT{mm};
        loca00  = PHI.local_MF{mm};
        ys1     = ys(loca0,:);
        Loca0 = PHI.Local_MF{mm};
        for ii = 1:d-e
            if kk == PHI.Omt
                theta1  = Theta(loca0,:)';
                if PHI.MATLABsolver
                    dys1    = MT_SHT_dy(data,ys1,theta1,Om1,repmat(p,1,PHI.pS2d),N);
                else
                    ys1(:,N/2+1:end)  = ys1(:,N/2+1:end)*Om1;
                    dys1    = MT_SHT_dy(data,ys1,theta1,Om1,repmat(p,1,PHI.pS2d),N);
                    dys1(:,N/2+1:end) = dys1(:,N/2+1:end)/Om1;
                end
                dYs1          = PHI.iGam*dys1;
                JY1(ii,Loca0) = 1e5 * reshape(dYs1,[],1);
            else
                Ys1     = Ys(loca00,:);
                JY1(ii,Loca0) = 1e5 * kron(eye(N),PHI.XI1_i{kk})*reshape(Ys1,[],1);
            end
            kk = kk+1;
        end
    end
end
%==========================================================
y = R1;
J = [JY1, Jom1, Jp1];
end