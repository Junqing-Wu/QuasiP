function [y,J] = MT_SHT_HT_QuasiP(data,PSI,ys,Om)
%MT_SHT_HT_QUASIP  Assemble hyper-time FSE-SHT boundary equations.
%
% This routine is the hyper-time counterpart of MT_SHT_MF_QuasiP. It matches
% the terminal state of each shooting segment to the initial state of the next
% segment and closes the trajectory family after one period of the selected
% base-frequency direction.
%
% For adjacent segments the residual is
%   ye_jj(ys_jj) - ys_jjp1 = 0.
% For d>1, the closing condition is
%   kron(eye(N),Gam*Omic*iGam)*ye_M - ys_1 = 0.
% Gam and iGam convert between trajectory samples and Fourier coefficients;
% Omic applies the phase shift caused by the remaining base frequencies.
%
% On input:
%   data : QuasiP data, including the hyper-time mesh and conversion matrices.
%   PSI  : Terminal states and optional sensitivities from the shooting solver.
%   ys   : Initial hyper-time trajectory samples for all segments.
%   Om   : Base-frequency vector; Omt selects the integration direction.
%
% On output:
%   y    : Hyper-time continuity and closure residual.
%   J    : Jacobian with respect to ys, Om, and p when PHI.isvar is true.
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
%----------------------------------------------------------
R0   = zeros(NL,1);         % NAE boundary residual R0
Jy0  = zeros(NL,NL);        % dR0/dys
Jom0 = zeros(NL,d);         % dR0/dOm
Jp0  = zeros(NL,a);         % dR0/dp
%==========================================================
if PHI.isvar
    ye     = PSI(:,1:N);
    ye_ys  = PSI(:,N+1:N+N^2);
    ye_Om  = PSI(:,N+N^2+1:N+N^2+N*d);
    ye_p   = PSI(:,N+N^2+N*d+1:end);
    %==========================================================
    %>>>>>>>>>>>>>>>>> yejj(ysjj) - ysjjp1 = 0 <<<<<<<<<<<<<<<<
    %==========================================================
    if d == 1
        N_jj = M_SHT;
    else
        N_jj = M_SHT-1;
    end
    for jj = 1:N_jj
        loca0 = PHI.local_HT{jj};
        Loca0 = PHI.Local_HT{jj};
        if jj ~= M_SHT
            loca1 = PHI.local_HT{jj+1};
            Loca1 = PHI.Local_HT{jj+1};
        else
            loca1 = PHI.local_HT{1};
            Loca1 = PHI.Local_HT{1};
        end
        % ---------------------------------------------
        ye_jj     = reshape(ye(loca0,:),[],1);
        ys_jjp1   = reshape(ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + ye_jj - ys_jjp1;
        % ---------------------------------------------
        ye_ys_jj  = MT_SHT_traj_yeys(ye_ys(loca0,:),PHI.pS2d,N);
        Jy0(Loca0,Loca0) = Jy0(Loca0,Loca0) + ye_ys_jj;
        Jy0(Loca0,Loca1) = Jy0(Loca0,Loca1) - eye(N*PHI.pS2d);
        % ---------------------------------------------
        ye_Om_jj        = reshape(ye_Om(loca0,:),[],d);
        for ii = 1:d
            Jom0(Loca0,ii) = Jom0(Loca0,ii) + ye_Om_jj(:,ii);
        end
        % ---------------------------------------------
        if data.isp
            ye_p_jj         = reshape(ye_p(loca0,:),[],a);
            Jp0(Loca0,:)    = Jp0(Loca0,:) + ye_p_jj;
        end
    end
    %==========================================================
    %>>>>>> (I_N \otimes G*Omic*iG)yeM(ysM) - ys1 = 0 <<<<<<<<<
    %==========================================================
    if d ~= 1
        [Omic,dOmic_dOm] = MT_SHT_Omic(PHI,Om1,Rho,d,Omt,Omnot);

        loca0 = PHI.local_HT{M_SHT};
        Loca0 = PHI.Local_HT{M_SHT};
        loca1 = PHI.local_HT{1};
        Loca1 = PHI.Local_HT{1};
        EGOG      = sparse( kron(eye(N),PHI.Gam*Omic*PHI.iGam) );
        % ---------------------------------------------
        ye_jj     = reshape(ye(loca0,:),[],1);
        ys_jjp1   = reshape(ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + EGOG*ye_jj - ys_jjp1;
        % ---------------------------------------------
        ye_ys_jj  = MT_SHT_traj_yeys(ye_ys(loca0,:),PHI.pS2d,N);
        Jy0(Loca0,Loca0) = Jy0(Loca0,Loca0) + EGOG*ye_ys_jj;
        Jy0(Loca0,Loca1) = Jy0(Loca0,Loca1) - eye(PHI.pS2d*N);
        % ---------------------------------------------
        ye_Om_jj        = reshape(ye_Om(loca0,:),[],d);
        for ii = 1:d
            EGO_omG        = sparse( kron(eye(N),PHI.Gam*dOmic_dOm{ii}*PHI.iGam) );
            Jom0(Loca0,ii) = EGOG*ye_Om_jj(:,ii) + EGO_omG*ye_jj;
        end
        % ---------------------------------------------
        if data.isp
            ye_p_jj         = reshape(ye_p(loca0,:),[],a);
            Jp0(Loca0,:)    = Jp0(Loca0,:) + EGOG*ye_p_jj;
        end
    end
    %==========================================================
    y = R0;
    J = [Jy0, Jom0, Jp0];
else
    ye     = PSI(:,1:N);
    %==========================================================
    %>>>>>>>>>>>>>>>>> yejj(ysjj) - ysjjp1 = 0 <<<<<<<<<<<<<<<<
    %==========================================================
    if d == 1
        N_jj = M_SHT;
    else
        N_jj = M_SHT-1;
    end
    for jj = 1:N_jj
        loca0 = PHI.local_HT{jj};
        Loca0 = PHI.Local_HT{jj};
        if jj ~= M_SHT
            loca1 = PHI.local_HT{jj+1};
        else
            loca1 = PHI.local_HT{1};
        end
        % ---------------------------------------------
        ye_jj     = reshape(ye(loca0,:),[],1);
        ys_jjp1   = reshape(ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + ye_jj - ys_jjp1;
    end
    %==========================================================
    %>>>>>> (I_N \otimes G*Omic*iG)yeM(ysM) - ys1 = 0 <<<<<<<<<
    %==========================================================
    if d ~= 1
        [Omic,~] = MT_SHT_Omic(PHI,Om1,Rho,d,Omt,Omnot);

        loca0 = PHI.local_HT{M_SHT};
        Loca0 = PHI.Local_HT{M_SHT};
        loca1 = PHI.local_HT{1};
        EGOG      = sparse( kron(eye(N),PHI.Gam*Omic*PHI.iGam) );
        % ---------------------------------------------
        ye_jj     = reshape(ye(loca0,:),[],1);
        ys_jjp1   = reshape(ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + EGOG*ye_jj - ys_jjp1;
    end
    %==========================================================
    y = R0;
    J = [];
end
end

% data1 = ye_jj(1:16);
% data2 = ys_jjp1(1:16);
% data3 = EGOG*ye_jj;
% data3 = data3(1:16);
% x = 1:16;
% width = 0.22;
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
%  xlabel('$\bar{\varphi}_2$','Interpreter','latex');
% ylabel('$y_1$','Interpreter','latex');
% leg = legend('$\bar{\mathbf{y}}_{e:3}$','$\bar{\mathbf{y}}_{s:1}$','$\hat{\mathbf{y}}_{e:3}$');
% leg.Interpreter = 'latex';
% grid on
% plot_setting_paper()
% xlim([0.5 16.5])
% xticks(1:16)
% hold off