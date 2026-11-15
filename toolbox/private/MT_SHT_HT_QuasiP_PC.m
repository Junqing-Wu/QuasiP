function [data, y, J] = MT_SHT_HT_QuasiP_PC(prob, data, u)
% HT-OPC relation:
%   OPC_i = <phase-direction derivative, continuation correction> = 0.
% The code evaluates this row using the initial hyper-time trajectory.
%MT_SHT_HT_QUASIP_PC  Build phase conditions for hyper-time shooting.
%
% Compute the orthogonality rows associated with the unknown hyper-time
% frequencies using the initial trajectory samples ys and the trajectory
% derivative MT_SHT_dy. These conditions remove the phase-shift nullspace
% before COCO continues the SHT solution.
%
% On input:
%   prob : COCO problem structure.
%   data : QuasiP model and PHI data.
%   u    : Shooting unknown vector.
%
% On output:
%   data : Updated callback data.
%   y    : Hyper-time shooting phase-condition residual.
%   J    : Its Jacobian with respect to the active unknowns.
PHI   = data.PHI;
NL    = PHI.nL;
N     = PHI.N;
d     = data.ode_opts.d;
e     = data.ode_opts.e;
a     = data.pdim;
%----------------------------------------------------------
R1   = zeros(d-e,1);        % OPC residual for unknown response frequencies
Jy1  = zeros(d-e,NL);       % dOPC/dC
Jom1 = zeros(d-e,d);        % dOPC/dOm
Jp1  = zeros(d-e,a);        % dOPC/dp
%----------------------------------------------------------
ys  = reshape(u(1:NL),[],N);
Om  = u(1+NL:d+NL);
p = u(d+NL+1:end);
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
%>>>>>>>> (I_N \otimes iG^\top*iG)dys1   orthogonal to ys1r <<<<<<<<<<
%>>>>>>>> (I_N \otimes iG^\top*Xi*iG)ys1 orthogonal to ys1r <<<<<<<<<<
%============================================================
if d > e
    for mm = 1:PHI.M
        kk = e + 1;
        loca0   = PHI.local_HT{mm};
        ys1     = ys(loca0,:);
        Loca0 = PHI.Local_HT{mm};
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
                Jy1(ii,Loca0) = 1e5* kron(eye(N),PHI.iGTiG)*reshape(dys1,[],1);
            else
                Jy1(ii,Loca0) = 1e5* kron(eye(N),PHI.iGTXI1iG{kk})*reshape(ys1,[],1);
            end
            kk = kk+1;
        end
    end
end
%==========================================================
y = R1;
J = [Jy1, Jom1, Jp1];
end