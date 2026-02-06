%% Core Equation Background (Quasi-Periodic Dynamics)
% 2ND ODEs: M(p)*X_,t,t(Om,t) +  F2(X(Om,t),X_,t(Om,t),om,t,p) = 0
% 1ST ODEs: B(p)*X_,t(Om,t)   =  F1(X(Om,t),om,t,p)
% 
% om = [Om1,Om2,...,Ome]^\top; (physical excitation frequencies)
% Om = [Om1,Om2,...,Ome,Omep1,...,Omd]^\top; (extended QP base frequencies, d >= e >= 0)
% 
% data.PHI.parametrization_tori = 'MP'  (Multi-periods parametrization for invariant tori)
% 
% MP Parametrization (Time → Phase Domain):
% 2ND ODEs: \sum_{i=1}^{d}\sum_{j=1}^{d} OmiOmj M(p)*X_,taui,tauj(Tau) + F2(X(Tau), \sum_{i=1}^{d}OmiX_,taui(Tau), tau, p ) = 0
% 1ST ODEs: \sum_{i=1}^{d} Omi B(p)*X_,taui(Tau) = F1(X(Tau), tau, p )
% 
% tau = [tau1,tau2,...,taue]^\top; (physical phase variables)
% Tau = [tau1,tau2,...,taue,tauep1,...,taud]^\top; (extended phase variables for MP tori)
% 
% Force Decomposition:
% 2ND ODEs: F2(...) = D(p)*X_,t + K(p)*X + N2(...) - E2(om,t,p) = 0
% 1ST ODEs: F1(...) = A(p)*X + N1(...) + E1(om,t,p)
% 
% 2ND → 1ST Order Conversion (Augmented State Z = [X; X_,t]):
% B = [D M; M 0]; A = [-K 0; 0 M]; N1 = [-N2(X); 0]; E1 = [E2; 0];
% OR B = I; A = [0 I; -iM*K -iM*D]; N = [0; -iM*N2(X)]; E = [0; iM*E2];

function [data, y, J] = MP_mVCF_QuasiP(prob, data, u)
    % MP_mVCF_QuasiP: Core solver for Multi-Periodic Variational Continuation Method (mVCF)
    %                 for quasi-periodic (QP) invariant tori in nonlinear dynamical systems
    % Inputs:
    %   - prob: Problem structure (QP dynamics settings, continuation parameters)
    %   - data: System data structure (DOF, PHI config, MP parametrization, force functions)
    %   - u: Optimization variable vector [Xd; Om; p]
    % Outputs:
    %   - data: Updated data structure (intermediate calculations, PHI config)
    %   - y: Combined residual vector [y1; y2] (main QP equilibrium + phase condition)
    %   - J: Jacobian matrix [J1_Xd J1_Om J1_p; J2_Xd J2_Om J2_p] (derivatives of y w.r.t u)

    %% 1. Extract Core System/Parametrization Parameters
    nL = data.nL;                  % Dimension of Xd (undetermined coefficients for QP tori)
    d = data.ode_opts.d;           % Number of extended QP base frequencies (d ≥ e)
    e = data.ode_opts.e;           % Number of physical excitation frequencies (e ≥ 0)
    n = data.DOF;                  % Degrees of freedom (DOF) of the dynamical system
    PHI = data.PHI;                % PHI config (MP torus parametrization, HB/CO/FD spectral methods)
    U = cell2mat(PHI.U);           % U-domain harmonic order counts (k_i)
    U_prod = prod(U);              % Total U-domain harmonic combinations (product of k_1^d)
    S = cell2mat(PHI.S);           % S-domain spectral component counts (bartau_i)
    S_prod = prod(S);              % Total S-domain spectral combinations (product of bartau_1^d)
    tau = PHI.tau';                % Transposed extended phase vector [τ₁, τ₂, ..., τ_d] (MP parametrization)

    %% 2. Decompose Optimization Variable Vector u
    Xd = u(1:nL);                  % Undetermined coefficients (slow-varying QP tori behavior)
    Xd_n = reshape(Xd,U_prod,n)';  % Reshape Xd to n×U_prod (DOF × U-domain harmonics)
    Om = u(nL+1:nL+d);             % Extended QP base frequency vector [Ω₁, Ω₂, ..., Ω_d]
    p  = u(nL+d+1:end);            % Physical parameter vector (stiffness/damping/mass)
    a = numel(p);                  % Dimension of physical parameter vector p

    %% 3. Initialize Residual/Jacobian for Constraints
    % 3.1 Main QP Equilibrium Constraint: Rd(Xd,Om,p) = 0 (core nonlinear system for fsolve)
    y1 = zeros(nL,1);               % Residual vector for Rd = 0 (dimension: nL)
    J1_Xd = zeros(nL,nL);           % Jacobian: ∂Rd/∂Xd (QP torus coefficients)
    J1_Om = zeros(nL,d);            % Jacobian: ∂Rd/∂Om (QP base frequencies)
    J1_p = zeros(nL,a);             % Jacobian: ∂Rd/∂p (physical parameters)

    % 3.2 Phase Condition Constraint: PC(Xd) = 0 (fix QP torus phase ambiguity)
    %     Dimension: d-e (gap between extended/physical frequencies)
    y2 = zeros(d-e,1);              % Residual vector for PC = 0
    J2_Xd = zeros(d-e,nL);          % Jacobian: ∂PC/∂Xd
    J2_Om = zeros(d-e,d);           % Jacobian: ∂PC/∂Om
    J2_p = zeros(d-e,a);            % Jacobian: ∂PC/∂p

    %% 4. Solve 1st/2nd Order QP ODEs
    switch data.order
        case{1} % 1st Order ODEs (State-Space Form for QP Dynamics)
            % 4.1 Compute Phase Derivative Matrices (Xi1 = ∑Ω_i*Xi.st{ii}, Xi0 = zero-order phase matrix)
            Xi1 = sparse(U_prod,U_prod);
            for ii = 1:d
                Xi1 = Xi1 + Om(ii)*PHI.Xi.st{ii};
            end
            Xi0 = PHI.Xi.zth;

            % 4.2 Map U-domain Xd to S-domain X0 (parallelized for efficiency)
            X0 = zeros(n,S_prod);
            parfor nn = 1:n
                X0(nn,:) = U2S_0(PHI,Xd_n(nn,:),U_prod);
            end

            % 4.3 Evaluate System Matrices/Forces (S-domain)
            B   = data.b(p);        % System matrix B(p) (1st-order ODE)
            B_P = data.bp(p);       % Derivative of B w.r.t physical parameters p
            F = data.f(tau,X0,p,'');% Nonlinear force F1(X,Tau,p) (S-domain)
            F_X = data.fx(tau,X0,p,'');% Derivative of F w.r.t X (S-domain)
            F_P = data.fp(tau,X0,p,'');% Derivative of F w.r.t p (S-domain)

            % 4.4 Compute Main Residual y1 (QP Equilibrium Constraint Rd = 0)
            %     y1 = -(B⊗Xi1)Xd + (A⊗Xi0)Xd + (I⊗Xi0)(Nd + Ed)
            BXi1 = kron(B,Xi1);
            y1 = y1 - BXi1*Xd;

            % Add linear term (A⊗Xi0)Xd
            if isfield(F, 'A')
                AXi0 = kron(F.A,Xi0);
                y1 = y1 + AXi0*Xd;
            end

            % Add nonlinear term Nd (S→U domain mapping)
            if isfield(F, 'N_tau')
                Nd_n = zeros(U_prod,n);
                non_zero_nn = find(~all(F.N_tau == 0, 2)); % Vectorized non-zero row filter
                for nn = non_zero_nn
                    Nd_n(:,nn) = S2U_1(PHI,F.N_tau(nn,:),S_prod);
                end
                Nd = reshape(Nd_n,[],1);
                if PHI.isCO1 == 0
                    y1 = y1 + Nd;
                else
                    IXi0 = kron(eye(n),Xi0);
                    y1 = y1 + IXi0*Nd;
                end
            end

            % Add excitation term Ed (S→U domain mapping)
            if isfield(F, 'E_tau')
                Ed_n = zeros(U_prod,n);
                non_zero_nn = find(~all(F.E_tau == 0, 2)); % Vectorized non-zero row filter
                for nn = non_zero_nn % Fixed: original "nn=1:non_zero_nn" (typo)
                    Ed_n(:,nn) = S2U_1(PHI,F.E_tau(nn,:),S_prod);
                end
                Ed = reshape(Ed_n,[],1);
                if PHI.isCO1 == 0
                    y1 = y1 + Ed;
                else
                    IXi0 = kron(eye(n),Xi0);
                    y1 = y1 + IXi0*Ed;
                end
            end

            % 4.5 Compute Jacobian J1_Om (∂Rd/∂Om)
            for ii = 1:d
                J1_Om(:,ii) = -kron(B,PHI.Xi.st{ii})*Xd;
            end

            % 4.6 Compute Jacobian J1_Xd (∂Rd/∂Xd)
            J1_Xd = J1_Xd - BXi1;
            if isfield(F_X, 'Ax_x')
                J1_Xd = J1_Xd + AXi0;
            end

            % Add nonlinear gradient Nd_Xd (S→U domain mapping)
            if isfield(F_X, 'N_x')
                Nd_Xd = zeros(nL,nL);
                is_nonzero = ~all(F_X.N_x == 0, 3); % 3D non-zero filter (ll,mm,S_prod)
                [non_zero_ll, non_zero_mm] = find(is_nonzero);
                for idx = 1:length(non_zero_ll)
                    ll = non_zero_ll(idx);
                    mm = non_zero_mm(idx);
                    Ndl_Xdm = S2U_2_withX(PHI, F_X.N_x(ll, mm, :), S_prod);
                    % Block-wise gradient accumulation (U_prod × U_prod per DOF pair)
                    row_range = (ll-1)*U_prod + 1 : ll*U_prod;
                    col_range = (mm-1)*U_prod + 1 : mm*U_prod;
                    Nd_Xd(row_range, col_range) = Nd_Xd(row_range, col_range) + Ndl_Xdm;
                end
                if PHI.isCO1 == 0
                    J1_Xd = J1_Xd + Nd_Xd;
                else
                    J1_Xd = J1_Xd + IXi0*Nd_Xd;
                end
            end

            % 4.7 Compute Jacobian J1_p (∂Rd/∂p)
            for ii = 1:a
                % Contribution from B_Pii (derivative of B w.r.t p(ii))
                B_Pii = B_P(:,:,ii);
                if ~isempty(find(B_Pii, 1)) % Efficient non-zero check (short-circuit)
                    J1_p(:,ii) = J1_p(:,ii) + kron(-B_Pii,Xi1)*Xd;
                end

                % Contribution from A_Pii (derivative of A w.r.t p(ii))
                if isfield(F_P, 'A_p')
                    A_Pii = F_P.A_p(:,:,ii);
                    if ~isempty(find(A_Pii, 1))
                        J1_p(:,ii) = J1_p(:,ii) + kron(A_Pii,Xi0)*Xd;
                    end
                end

                % Contribution from NE_p (nonlinear+excitation derivative w.r.t p(ii))
                if isfield(F_P,'NE_p')
                    NE_p_slice = F_P.NE_p(:,ii,:);
                    NE0_Pii = reshape(NE_p_slice,n,S_prod);
                    non_zero_nn = find(~all(NE0_Pii == 0, 2));
                    if ~isempty(non_zero_nn)
                        NEd_Pii = zeros(U_prod,n);
                        for nn = non_zero_nn % Fixed: original "nn=1:non_zero_nn" (typo)
                            NEd_Pii(:,nn) = S2U_1(PHI,NE0_Pii(nn,:),S_prod);
                        end
                        NEd_Pii = reshape(NEd_Pii,[],1);
                        if PHI.isCO1 == 0
                            J1_p(:,ii) = J1_p(:,ii) + NEd_Pii;
                        else
                            J1_p(:,ii) = J1_p(:,ii) + IXi0*NEd_Pii;
                        end
                    end
                end
            end

            % 4.8 Compute Phase Condition Jacobian J2_Xd (fix QP torus phase ambiguity)
            if d>e
               kk = e+1;
               for ii = 1:d-e
                   ILambda_ii = kron(eye(n),PHI.Lambda{1,kk});
                   J2_Xd(ii,:) = transpose(ILambda_ii*Xd);
                   kk = kk+1;
               end
            end

            % 4.9 Assemble Total Residual/Jacobian
            y = [y1;y2];
            J = [J1_Xd,J1_Om,J1_p
                 J2_Xd,J2_Om,J2_p];

        case{2} % 2nd Order ODEs (To be implemented)
    end

end

%% Helper Function 1: U2S_0 - Map U-domain Xd to S-domain X0 (spectral projection)
function X0 = U2S_0(PHI,Xd,U_prod)
    % U2S_0: U-domain (harmonic) → S-domain (spectral) mapping for QP dynamics
    % Inputs:
    %   - PHI: MP parametrization config (HB/CO/FD methods, U/S dimensions)
    %   - Xd: U-domain coefficient vector (dimension: U_prod)
    %   - U_prod: Total U-domain harmonic combinations
    % Output:
    %   - X0: S-domain spectral vector (dimension: S_prod)

    Siip1Td = 1;
    U1Tii = U_prod;
    U = cell2mat(PHI.U);
    S = cell2mat(PHI.S);

    Xii = Xd;
    % Iterative spectral projection (backward over extended frequencies)
    for kk = PHI.d:-1:1
        U1Tii = U1Tii/U(kk);
        % Reshape + transpose for spectral transformation
        Xii_wave = pagetranspose( reshape(Xii,Siip1Td,U(kk),U1Tii) );
        Xii_wave = reshape(Xii_wave,U(kk),Siip1Td*U1Tii);
        
        % Spectral method (HB/CO/FD)
        switch PHI.method{kk}
            case{'HB'} % Harmonic Balance (HB)
                GAM = PHI.Gamma_i{1,kk};
                Xiiminus1_bar = GAM*Xii_wave;
            case{'CO'} % Collocation (CO)
                Xiiminus1_bar = Xii_wave;
            case{'FD'} % Finite Difference (FD)
                Xiiminus1_bar = Xii_wave;
        end
        
        Siip1Td = Siip1Td*S(kk);
        Xii = reshape( Xiiminus1_bar, Siip1Td, U1Tii  );
    end
    X0 = Xii';
end

%% Helper Function 2: S2U_1 - Map S-domain N0 to U-domain Nd (inverse spectral projection)
function Nd = S2U_1(PHI,N0,S_prod)
    % S2U_1: S-domain (spectral) → U-domain (harmonic) mapping for nonlinear terms
    % Inputs:
    %   - PHI: MP parametrization config (HB/CO/FD methods, U/S dimensions)
    %   - N0: S-domain nonlinear force vector (dimension: S_prod)
    %   - S_prod: Total S-domain spectral combinations
    % Output:
    %   - Nd: U-domain nonlinear force matrix (dimension: U_prod × n)

    U = cell2mat(PHI.U);
    S = cell2mat(PHI.S);

    S_iiTd = S_prod;
    U_1Tiiminus1 = 1;
    Niiminus1 = N0;

    % Iterative inverse spectral projection (forward over extended frequencies)
    for kk = 1:PHI.d
        S_iiTd = S_iiTd/S(kk);
        % Reshape + transpose for inverse spectral transformation
        Niiminus1_bar = pagetranspose(reshape(Niiminus1,U_1Tiiminus1,S(kk),S_iiTd) );
        Niiminus1_bar = reshape(Niiminus1_bar,S(kk),U_1Tiiminus1*S_iiTd);

        % Inverse spectral method (HB/CO/FD)
        switch PHI.method{kk}
            case{'HB'} % Inverse Harmonic Balance
                iGAM = PHI.Gamma_i{3,kk};
                Fnlii = iGAM*Niiminus1_bar;
            case{'CO'} % Collocation (identity mapping)
                Fnlii = Niiminus1_bar;
            case{'FD'} % Finite Difference (identity mapping)
                Fnlii = Niiminus1_bar;
        end
        
        U_1Tiiminus1 = U_1Tiiminus1*U(kk);
        Niiminus1 = reshape(Fnlii,U_1Tiiminus1,S_iiTd);
    end
    Nd = Niiminus1;
end

%% Helper Function 3: S2U_2_withX - Map S-domain N0_X0 to U-domain Nd_Xd (gradient projection)
function Nd_Xd = S2U_2_withX(PHI,N0_X0,S_prod)
    % S2U_2_withX: S-domain → U-domain mapping for nonlinear force gradients (∂N/∂Xd)
    % Inputs:
    %   - PHI: MP parametrization config (HB/CO/FD methods, U/S dimensions)
    %   - N0_X0: S-domain gradient tensor (dimension: 1×S_prod)
    %   - S_prod: Total S-domain spectral combinations
    % Output:
    %   - Nd_Xd: U-domain gradient matrix (dimension: U_prod × U_prod)

    U = cell2mat(PHI.U);
    S = cell2mat(PHI.S);
    S_iiTd = S_prod;
    U_1Tiiminus1 = 1;

    N_Xiiminus1 = reshape(N0_X0,1,[]);

    % Iterative gradient projection (forward over extended frequencies)
    for kk = 1:PHI.d
        S_iiTd = S_iiTd/S(kk);
        % Reshape + transpose for gradient transformation
        N_Xiiminus1_bar = pagetranspose( reshape(N_Xiiminus1,U_1Tiiminus1^2,S(kk),S_iiTd) );
        N_Xiiminus1_bar = reshape( N_Xiiminus1_bar, S(kk),1,U_1Tiiminus1^2*S_iiTd );

        % Gradient transformation (HB/CO/FD)
        switch PHI.method{kk}
            case{'HB'} % Harmonic Balance gradient
                GAM = PHI.Gamma_i{1,kk};
                iGAM = PHI.Gamma_i{3,kk};
                dFnl_dZiiminus1_bar_SU = repmat(N_Xiiminus1_bar,1,U(kk));
                dFnl_dZii = dFnl_dZiiminus1_bar_SU.*GAM;
                dFnl_dZii = pagemtimes( iGAM,dFnl_dZii );
            case{'CO'} % Collocation gradient (diagonal mapping)
                dFnl_dZiiminus1_bar_SS = zeros(S(kk),S(kk),U_1Tiiminus1^2*S_iiTd);
                for jj = 1:U_1Tiiminus1^2*S_iiTd
                    dFnl_dZiiminus1_bar_SS(:,:,jj) = diag(N_Xiiminus1_bar(:,:,jj));
                end
                dFnl_dZii = dFnl_dZiiminus1_bar_SS;
            case{'FD'} % Finite Difference gradient (diagonal mapping)
                dFnl_dZiiminus1_bar_SS = zeros(S(kk),S(kk),U_1Tiiminus1^2*S_iiTd);
                for jj = 1:U_1Tiiminus1^2*S_iiTd
                    dFnl_dZiiminus1_bar_SS(:,:,jj) = diag(N_Xiiminus1_bar(:,:,jj));
                end
                dFnl_dZii = dFnl_dZiiminus1_bar_SS;
        end
        
        % Reshape + transpose for next iteration
        N_Xiiminus1 = pagetranspose( reshape(dFnl_dZii,U(kk),U(kk)*U_1Tiiminus1,U_1Tiiminus1,S_iiTd) );
        U_1Tiiminus1 = U_1Tiiminus1*U(kk);
        N_Xiiminus1 = reshape(N_Xiiminus1,U_1Tiiminus1^2,S_iiTd);
    end
    Nd_Xd = transpose( reshape(N_Xiiminus1,U_1Tiiminus1,U_1Tiiminus1) );
end