function [y,J] = MT_SHT_MF_QuasiP(data,PSI,Ys,Om)
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
R0   = zeros(NL,1);         % 非线性代数方程组Rd的残差
JY0  = zeros(NL,NL);        % Rd关于Xd的导数
Jom0 = zeros(NL,d);         % Rd关于Om的导数
Jp0  = zeros(NL,a);         % Rd关于p 的单数
%==========================================================
if PHI.isvar
    ye     = PSI(:,1:N);
    ye_ys  = PSI(:,N+1:N+N^2);
    ye_Om  = PSI(:,N+N^2+1:N+N^2+N*d);
    ye_p   = PSI(:,N+N^2+N*d+1:end);
    %==========================================================
    %>>>>>>>>>>>>>>>>> Yejj(ysjj) - Ysjjp1 = 0 <<<<<<<<<<<<<<<<
    %==========================================================
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
    %==========================================================
    %>>>>>>>>> (I_N \otimes Omic)YeM(YsM) - Ys1 = 0 <<<<<<<<<<<
    %==========================================================
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
    %==========================================================
    y = R0;
    J = [JY0, Jom0, Jp0];
else
    ye     = PSI(:,1:N);
    %==========================================================
    %>>>>>>>>>>>>>>>>> Yejj(ysjj) - Ysjjp1 = 0 <<<<<<<<<<<<<<<<
    %==========================================================
    for jj = 1:M_SHT-1
        loca0 = PHI.local_HT{jj};
        loca1 = PHI.local_MF{jj+1};
        Loca0 = PHI.Local_MF{jj};
        % ---------------------------------------------
        Ye_jj     = reshape(PHI.iGam*ye(loca0,:),[],1);
        Ys_jjp1   = reshape(Ys(loca1,:),[],1);
        R0(Loca0) = R0(Loca0) + Ye_jj - Ys_jjp1;
    end
    %==========================================================
    %>>>>>>>>> (I_N \otimes Omic)YeM(YsM) - Ys1 = 0 <<<<<<<<<<<
    %==========================================================
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
    %==========================================================
    y = R0;
    J = [];
end
end