function [y,J] = MT_SHT_MF_QuasiP(data,PSI,Ys,Om)
%MT_SHT_MF_QUASIP  Assemble the multi-frequency FSE-SHT boundary equations.
%
% This function is the boundary-condition part of the multi-trajectory,
% Fourier-series expansion-based shooting formulation. MT_SHT_QuasiP first
% integrates every trajectory segment and stores the terminal states and
% sensitivities in PSI. This function then converts those terminal states to
% the multi-frequency coefficient domain and imposes the continuity equations.
%
% For adjacent trajectory segments, the required relation is
%   Ye_jj(ys_jj) - Ys_jjp1 = 0,   jj = 1,...,M_SHT-1.
% For the final segment, the initial section is reached after a phase shift
% caused by the remaining base frequencies. The closure relation is
%   kron(eye(N),Omic)*Ye_M(ys_M) - Ys_1 = 0.
% Here Omic is constructed by MT_SHT_Omic and contains the period factors in
% PHI.P, the frequency ratios Rho, and the Fourier index combinations PHI.kQ.
%
% The residual vector returned as y is ordered according to PHI.Local_MF.
% When PHI.isvar is true, J contains the derivatives of the same boundary
% equations with respect to the initial coefficients Ys, the base frequencies
% Om, and the model parameters p.
%
% [Y,J] = MT_SHT_MF_QUASIP(DATA, PSI, YS, OM)
%
% On input:
%
% DATA : QuasiP data structure. DATA.PHI supplies the segment maps, the
%        Fourier conversion matrices PHI.Gam and PHI.iGam, the period factors,
%        and the selected shooting mesh. DATA.isp indicates whether p is active.
%
% PSI  : Terminal trajectory data. Its first N columns are terminal states.
%        If PHI.isvar is true, the following columns store derivatives with
%        respect to Ys, Om, and p in the order used by MT_SHT_QuasiP.
%
% YS   : Initial multi-frequency coefficient blocks for the M_SHT trajectory
%        segments. The block locations are given by PHI.local_MF and
%        PHI.Local_MF.
%
% OM   : Current vector of d base frequencies. Omt selects the integration
%        direction, Om1=Om(Omt), and Rho=Om(Omnot)/Om1 gives the relative
%        drift of the remaining frequencies during one integration period.
%
% On output:
%
% Y    : Boundary residual vector. A converged shooting torus satisfies y=0.
%
% J    : Jacobian of y with respect to the active unknowns. It is empty when
%        PHI.isvar is false.
%
% See also: MT_SHT_QuasiP, MT_SHT_HT_QuasiP, MT_SHT_Omic,
%           MT_SHT_MF_QuasiP_PC.
PHI   = data.PHI;
NL    = PHI.nL;
N     = PHI.N;
d     = data.ode_opts.d;
a     = data.pdim;
M_SHT = PHI.M;
Omt   = PHI.Omt;
Omnot = PHI.Omnot;
Om1   = Om(Omt);
Om2d  = Om(Omnot);
Rho   = Om2d./Om1;
% Basic dimensions and frequency ratios used by all segment equations.
R0   = zeros(NL,1);         % Boundary residual vector R0
JY0  = zeros(NL,NL);        % dR0/dYs
Jom0 = zeros(NL,d);         % dR0/dOm
Jp0  = zeros(NL,a);         % dR0/dp

% Allocate the residual and Jacobian blocks. The three Jacobian blocks are
% d y/d Ys, d y/d Om, and d y/d p.
%==========================================================
if PHI.isvar
    ye     = PSI(:,1:N);
    ye_ys  = PSI(:,N+1:N+N^2);
    ye_Om  = PSI(:,N+N^2+1:N+N^2+N*d);
    ye_p   = PSI(:,N+N^2+N*d+1:end);
    % Internal segment matching: the terminal coefficients of segment jj must equal
    % the initial coefficients of segment jj+1.
    for jj = 1:M_SHT-1
        loca0 = PHI.local_HT{jj};
        loca1 = PHI.local_MF{jj+1};

        Loca0 = PHI.Local_MF{jj};
        Loca1 = PHI.Local_MF{jj+1};
        % ---------------------------------------------
        Ye_jj     = reshape(PHI.iGam*ye(loca0,:),[],1);
        Ys_jjp1   = reshape(Ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + Ye_jj - Ys_jjp1;
        % ---------------------------------------------
        ye_ys_jj  = MT_SHT_traj_yeys(ye_ys(loca0,:),PHI.pS2d,N);
        JY0(Loca0,Loca0) = JY0(Loca0,Loca0) + kron(eye(N),PHI.iGam)*ye_ys_jj*kron(eye(N),PHI.Gam);
        JY0(Loca0,Loca1) = JY0(Loca0,Loca1) - eye(N*PHI.U);
        % ---------------------------------------------
        Ye_Om_jj        = reshape(PHI.iGam*ye_Om(loca0,:),[],d);
        for ii = 1:d
            Jom0(Loca0,ii) = Jom0(Loca0,ii) + Ye_Om_jj(:,ii);
        end
        % ---------------------------------------------
        if data.isp
            Ye_p_jj         = reshape(PHI.iGam*ye_p(loca0,:),[],a);
            Jp0(Loca0,:)    = Jp0(Loca0,:) + Ye_p_jj;
        end
    end
    % Final torus closure: apply the Fourier phase-shift matrix to the terminal
    % coefficients of the last segment and match them to the first segment.
    if d ~= 1
        [Omic,dOmic_dOm] = MT_SHT_Omic(PHI,Om1,Rho,d,Omt,Omnot);

        loca0 = PHI.local_HT{M_SHT};
        loca1 = PHI.local_MF{1};
        Loca0 = PHI.Local_MF{M_SHT};
        Loca1 = PHI.Local_MF{1};
        EO    = sparse( kron(eye(N),Omic) );
        % ---------------------------------------------
        Ye_jj     = reshape(PHI.iGam*ye(loca0,:),[],1);
        Ys_jjp1   = reshape(Ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + EO*Ye_jj - Ys_jjp1;
        % ---------------------------------------------
        ye_ys_jj  = MT_SHT_traj_yeys(ye_ys(loca0,:),PHI.pS2d,N);
        JY0(Loca0,Loca0) = JY0(Loca0,Loca0) + EO*PHI.EiG*ye_ys_jj*PHI.EG;
        JY0(Loca0,Loca1) = JY0(Loca0,Loca1) - eye(N*PHI.U);
        % ---------------------------------------------
        Ye_Om_jj        = reshape(PHI.iGam*ye_Om(loca0,:),[],d);
        for ii = 1:d
            EO_om          = sparse( kron(eye(N),dOmic_dOm{ii}) );
            Jom0(Loca0,ii) = EO*Ye_Om_jj(:,ii) + EO_om*Ye_jj;
        end
        % ---------------------------------------------
        if data.isp
            Ye_p_jj         = reshape(PHI.iGam*ye_p(loca0,:),[],a);
            Jp0(Loca0,:)    = Jp0(Loca0,:) + EO*Ye_p_jj;
        end
    end
    y = R0;
    J = [JY0, Jom0, Jp0];
else
    ye     = PSI(:,1:N);
    % Internal segment matching: the terminal coefficients of segment jj must equal
    % the initial coefficients of segment jj+1.
    for jj = 1:M_SHT-1
        loca0 = PHI.local_HT{jj};
        loca1 = PHI.local_MF{jj+1};
        Loca0 = PHI.Local_MF{jj};
        % ---------------------------------------------
        Ye_jj     = reshape(PHI.iGam*ye(loca0,:),[],1);
        Ys_jjp1   = reshape(Ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + Ye_jj - Ys_jjp1;
    end
    % Final torus closure: apply the Fourier phase-shift matrix to the terminal
    % coefficients of the last segment and match them to the first segment.
    if d ~= 1
        [Omic,~] = MT_SHT_Omic(PHI,Om1,Rho,d,Omt,Omnot);
        loca0 = PHI.local_HT{M_SHT};
        loca1 = PHI.local_MF{1};
        Loca0 = PHI.Local_MF{M_SHT};
        EO    = sparse( kron(eye(N),Omic) );
        % ---------------------------------------------
        Ye_jj     = reshape(PHI.iGam*ye(loca0,:),[],1);
        Ys_jjp1   = reshape(Ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + EO*Ye_jj - Ys_jjp1;
    end
    y = R0;
    J = [];
end
end


% data1 = Ye_jj(1:7);
% data2 = Ys_jjp1(1:7);
% data3 = EO*Ye_jj;
% data3 = data3(1:7);
% x = 1:7;
% width = 0.25;
% figure('Color','w');
% b1 = bar(x-width, data1, width);
% b1.FaceColor = [0 0 1];
% b1.FaceAlpha = 0.5;
% hold on
% b2 = bar(x, data2, width);
% b2.FaceColor = [1 0 0];
% b2.FaceAlpha = 0.5;
% b3 = bar(x+width, data3, width);
% b3.FaceColor = [0.6 0 0.8];
% b3.FaceAlpha = 0.5;
%  xlabel('$k_2$','Interpreter','latex');
% ylabel('$y_1$','Interpreter','latex');
% leg = legend('${\mathbf{Y}}_{e:3}$','${\mathbf{Y}}_{s:1}$','$\hat{\mathbf{Y}}_{e:3}$');
% leg.Interpreter = 'latex';
% grid on
% plot_setting_paper()
% xlim([0.5 7.5])
% ylim([-6 2])
% xticks(1:7)
% hold off