function [y,J] = MT_SHT_HT_QuasiP(data,PSI,ys,Om)
PHI   = data.PHI;
NL    = PHI.nL;
N     = PHI.N;
d     = data.ode_opts.d;
a     = data.pdim;
M_SHT = PHI.M_SHT;
Omt   = PHI.Omt;
Omnot = PHI.Omnot;
Om1   = Om(Omt);
Om2d  = Om(Omnot);
Rho   = Om2d./Om1;
%----------------------------------------------------------
R0   = zeros(NL,1);         % 非线性代数方程组R的残差
Jy0  = zeros(NL,NL);        % R关于y0的导数
Jom0 = zeros(NL,d);         % R关于Om的导数
Jp0  = zeros(NL,a);         % R关于p 的单数
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