function PHI = MP_mVCF_COM_PHI(PHI_INITIAL)
% MP_mVCF_COM_PHI - Core operator construction for Multi-step Variable Coefficient Formula (m-VCF)
%                   Integrates Alternating U-S Domain Method (AUS) and Phase Condition (PC) logic
%                   Supports HB (Harmonic Balance)/FD (Finite Difference)/CO (Collocation) methods
%                   Critical Note: CO with Gauss points requires special adjustment of invariant matrices
%
% INPUT:
%   PHI_INITIAL - Initial configuration structure for m-VCF:
%                 .d          = Number of steps (dimensions) in m-VCF
%                 .U{ii}      = Number of discrete periodic points in U-domain for step ii
%                 .k{ii}      = Method-specific parameter for step ii (CO: [poly_degree, point_type(1=Gauss)])
%                 .S{ii}      = Number of S-domain sampling points for Alternating U-S Domain Method (AUS) (step ii)
%                 .S_wave{ii} = Number of high-resolution sampling points for post-processing refinement (step ii)
%                 .method{ii} = Numerical method for step ii: 'HB'/'FD'/'CO'
%
% OUTPUT:
%   PHI - Updated structure containing core m-VCF operators:
%         .Xi_i      = Step-wise invariant matrices (per-step sub-matrices, core of m-VCF)
%         .Xi        = Global combined matrix (joint concatenation of invariant matrices across all steps)
%         .tau       = Multi-dimensional periodic grid points (basis for U↔S domain mapping)
%         .Lambda    = Weighted operator matrices dedicated to Phase Condition (PC)
%         .Gamma_i   = Interpolation matrices dedicated to Alternating U-S Domain Method (AUS) (U↔S mapping)
%         .W_Gamma_i = High-resolution interpolation matrices dedicated to post-processing refinement

%% 1. Initialize m-VCF structure and extract core parameters
PHI = PHI_INITIAL;          % Inherit initial m-VCF configuration
d = PHI_INITIAL.d;          % Number of steps (dimensions) in m-VCF

%% 2. Step-wise invariant matrices (core sub-matrices for m-VCF)
% Xi_i: Step-wise invariant matrices (3×d cell array) - unadjusted per-step core operators
%   Xi_i{1,ii} = 0th-order invariant matrix for step ii (dimension: U_i×U_i)
%   Xi_i{2,ii} = 1st-order invariant matrix for step ii (dimension: U_i×U_i)
%   Xi_i{3,ii} = 2nd-order invariant matrix for step ii (dimension: U_i×U_i)
%   Core unmodified operators for each step in m-VCF
Xi_i      = cell(3,d);     

% W_Xi_i: Adjusted invariant matrices for CO (Gauss points) (2×d cell array)
%   W_Xi_i{1,ii} = Adjusted 1st-order matrix for step ii (dimension: U_i×U_i)
%   W_Xi_i{2,ii} = Adjusted 2nd-order matrix for step ii (dimension: U_i×U_i)
%   Critical adjustment: Required for CO with Gauss points (mandatory for AUS/PC implementation)
%   Note: W_Xi_i ≠ Xi_i for CO (Gauss points), identical for HB/FD
W_Xi_i    = cell(2,d);

%% 3. Global combined matrix (integrated across all m-VCF steps)
% Xi: Global combined matrix - aggregated operators from all step-wise invariant matrices
%   Xi.zth = Multi-dimensional 0th-order global matrix (dimension: U_1^d × U_1^d)
%   Xi.st  = 1st-order global matrix per step (1×d cell array, dimension: U_1^d × U_1^d)
%   Xi.nd  = 2nd-order global matrix per step pair (d×d cell array, dimension: U_1^d × U_1^d)
Xi.zth    = 1;            % Initialize 0th-order global matrix (identity seed for kron product)
Xi.st     = cell(1,d);    % 1st-order global matrix (step-specific differentiation)
Xi.nd     = cell(d,d);    % 2nd-order global matrix (step-pair specific differentiation)

%% 4. Specialized operator matrices for PC/AUS/post-processing
% Lambda: Weighted multi-dimensional operators for Phase Condition (PC) (2×d cell array)
%   Lambda{1,ii} = Weighted 1st-order operator for PC (dimension: nU_1^d × nU_1^d)
%   Lambda{2,ii} = Weighted 2nd-order operator for PC (dimension: nU_1^d × nU_1^d)
Lambda    = cell(2,d);   

% Gamma_i: Interpolation matrices for Alternating U-S Domain Method (AUS) (3×d cell array)
%   Gamma_i{1,ii} = 0th-order interpolation matrix (U→S mapping, dimension: S_1^d × U_1^d)
%   Gamma_i{2,ii} = 1st-order derivative interpolation matrix (AUS, dimension: S_1^d × U_1^d)
%   Gamma_i{3,ii} = Inverse 0th-order interpolation matrix (S→U mapping, dimension: U_1^d × S_1^d)
Gamma_i   = cell(3,d);   

% W_Gamma_i: High-resolution interpolation matrices for post-processing refinement (2×d cell array)
%   W_Gamma_i{1,ii} = High-res 0th-order interpolation matrix (dimension: Sw_1^d × U_1^d)
%   W_Gamma_i{2,ii} = High-res 1st-order derivative interpolation matrix (dimension: Sw_1^d × U_1^d)
W_Gamma_i   = cell(2,d);   

%% 5. Extract total discrete points across all steps (U-domain)
U = cell2mat(PHI.U);  % U = [U_1, U_2, ..., U_d] - discrete points per step (U-domain)

%% 6. Construct step-wise operators (loop over each m-VCF step)
for ii = 1:d
    % Extract step-specific parameters
    Ui = PHI.U{ii};         % U-domain discrete points for step ii
    ki = PHI.k{ii};         % Method parameter for step ii (CO: [poly_degree, flag])
    Si = PHI.S{ii};         % S-domain sampling points for AUS (step ii)
    Swi = PHI.S_wave{ii};   % High-res sampling points for post-processing (step ii)
    
    % Construct operators based on selected numerical method
    switch PHI.method{ii}
        case{'HB'}
            % Harmonic Balance (HB) - unadjusted invariant matrices
            [Xi_i{2,ii},Xi_i{3,ii},Gamma_i{1,ii},Gamma_i{3,ii},W_Gamma_i{1,ii}] = HB(ki,Ui,Si,Swi);
            Xi_i{1,ii} = eye(Ui);                  % 0th-order HB operator (identity matrix)
            Gamma_i{2,ii} = Gamma_i{1,ii}*Xi_i{2,ii};  % 1st-order AUS interpolation (HB)
            W_Gamma_i{2,ii} = W_Gamma_i{1,ii}*Xi_i{2,ii};  % High-res 1st-order interpolation (HB)
            W_Xi_i(1:2,ii) = Xi_i(2:3,ii);         % W_Xi_i = Xi_i for HB (no adjustment)
            
        case{'FD'}
            % Finite Difference (FD) - unadjusted invariant matrices
            [Xi_i{2,ii},Xi_i{3,ii},W_Gamma_i{1,ii},W_Gamma_i{2,ii}] = FD(ki,Ui,Si,Swi);
            Xi_i{1,ii} = eye(Ui);                  % 0th-order FD operator (identity matrix)
            Gamma_i{1,ii} = eye(Ui);               % 0th-order AUS interpolation (FD, identity)
            Gamma_i{2,ii} = Xi_i{2,ii};            % 1st-order AUS derivative (FD)
            Gamma_i{3,ii} = eye(Ui);               % Inverse AUS interpolation (FD, identity)
            W_Xi_i(1:2,ii) = Xi_i(2:3,ii);         % W_Xi_i = Xi_i for FD (no adjustment)
            
        case{'CO'}
            % Collocation (CO) - adjusted invariant matrices (critical for Gauss points)
            [Xi_i(1:3,ii),W_Xi_i(1:2,ii),W_Gamma_i{1,ii},W_Gamma_i{2,ii}] = CO(ki,Ui,Si,Swi);
            Gamma_i{1,ii} = eye(Ui);               % 0th-order AUS interpolation (CO, identity)
            Gamma_i{2,ii} = W_Xi_i{1,ii};          % 1st-order AUS derivative (CO, adjusted)
            Gamma_i{3,ii} = eye(Ui);               % Inverse AUS interpolation (CO, identity)
    end
    
    % Update global 0th-order matrix (kron product across steps)
    Xi.zth = sparse( kron(Xi.zth,Xi_i{1,ii}) );

    % Calculate dimension separation factors for kronecker product
    U1_im1 = prod(U(1:ii-1));  % Product of U_1 to U_{ii-1} (preceding steps)
    Uip1_d = prod(U(ii+1:d));  % Product of U_{ii+1} to U_d (succeeding steps)
    
    % Construct PC-specific weighted global operators (Lambda)
    Lambda{1,ii} = sparse( kron( kron(eye(U1_im1),W_Xi_i{1,ii}), eye(Uip1_d)) );
    Lambda{2,ii} = sparse( kron( kron(eye(U1_im1),W_Xi_i{2,ii}), eye(Uip1_d)) );
end

%% 7. Construct global 1st-order operators (Xi.st)
% Xi.st{ii}: 1st-order differentiation along step ii (0th-order for other steps)
for ii =1:d
    Xi.st{ii} = 1;  % Initialize kron product accumulator
    for jj = 1:d
        if ii == jj
            % Use 1st-order invariant matrix for target step ii
            Xi.st{ii} = sparse( kron(Xi.st{ii},Xi_i{2,jj}) );
        else
            % Use 0th-order invariant matrix for non-target steps
            Xi.st{ii} = sparse( kron(Xi.st{ii},Xi_i{1,jj}) );            
        end
    end
end

%% 8. Construct global 2nd-order operators (Xi.nd)
% Xi.nd{ii,jj}: 2nd-order differentiation (pure ii or mixed ii/jj)
for ii = 1:d
    for jj = 1:d
        Xi.nd{ii,jj} = 1;  % Initialize kron product accumulator
        for kk = 1:d
            if kk == ii || kk == jj
                if ii == jj
                    % Pure 2nd-order differentiation (same step: use 2nd-order matrix)
                   Xi.nd{ii,jj} = sparse( kron(Xi.nd{ii,jj},Xi_i{3,kk}) ); 
                else
                    % Mixed 2nd-order differentiation (different steps: use 1st-order matrix)
                   Xi.nd{ii,jj} = sparse( kron(Xi.nd{ii,jj},Xi_i{2,kk}) ); 
                end
            else
                % Use 0th-order matrix for non-target steps
                Xi.nd{ii,jj} = sparse( kron(Xi.nd{ii,jj},Xi_i{1,kk}) ); 
            end
        end
    end
end

%% 9. Generate multi-dimensional periodic grid points (tau)
% tau: [S_1^d × d] matrix - U/S domain grid points (each row = (τ₁, τ₂, ..., τ_d))
S   = cell2mat (PHI.S);  % S = [S_1, S_2, ..., S_d] - AUS sampling points per step (S-domain)
% Initialize 1D grid for first step (U/S domain mapping basis)
tau_1im1 = (0:2*pi/S(1):2*pi-2*pi/S(1))';
S_1im1 = S(1);

% Extend to multi-dimensional grid via kronecker product (U↔S mapping)
for ii = 2:d
    tau_i = (0:2*pi/S(ii):2*pi-2*pi/S(ii))';  % 1D grid for step ii (S-domain)
    % Combine existing grid with new step (preserve periodicity)
    tau_1im1 = [ kron( ones(S(ii),1),tau_1im1 ) , kron( tau_i,ones(S_1im1,1) ) ];
    S_1im1 = S_1im1*S(ii);  % Update total AUS sampling points (S-domain)
end
tau = tau_1im1;  % Final multi-dimensional periodic grid (U↔S domain)

%% 10. Assign all operators to output structure
PHI.Xi_i = Xi_i;          % Step-wise invariant matrices (m-VCF core)
PHI.Xi = Xi;              % Global combined matrix (integrated operators)
PHI.tau = tau;            % Multi-dimensional periodic grid points
PHI.Lambda = Lambda;      % PC-dedicated weighted operators
PHI.Gamma_i = Gamma_i;    % AUS-dedicated interpolation matrices (U↔S)
PHI.W_Gamma_i = W_Gamma_i;% Post-processing refinement matrices

PHI.L = prod(cell2mat(PHI.U));
end



function [Xi1,Xi2,H_iDFT,H_DFT,W_H_iDFT] = HB(k_HB,L_HB,S_HB,SW_HB)
% HB - Core implementation of Harmonic Balance (HB) method
%
% INPUTS:
%   k_HB    - Harmonic order (frequency multiplier for HB method)
%   L_HB    - Total number of harmonic components (1 + 2*number_of_harmonics)
%   S_HB    - Number of sampling points for discrete Fourier transform (DFT)
%   SW_HB   - Number of sampling points for weighted DFT (for post-processing/verification)
%
% OUTPUTS:
%   Xi1     - First-order harmonic operator matrix (skew-symmetric diagonal)
%   Xi2     - Second-order harmonic operator matrix (square of Xi1)
%   H_iDFT  - Inverse DFT matrix for HB sampling (frequency domain -> time domain)
%   H_DFT   - DFT matrix for HB sampling (time domain -> frequency domain)
%   W_H_iDFT- Weighted inverse DFT matrix (higher-resolution sampling)
%
% DESCRIPTION:
%   This function constructs core matrices for Harmonic Balance method, including
%   harmonic operator matrices (Xi1/Xi2) and DFT/inverse DFT matrices for
%   periodic solution reconstruction.

%% 1. Construct harmonic operator vectors
% v1/v2: Auxiliary vectors for building skew-symmetric harmonic operator matrices
v1 = zeros(1,L_HB-1);          % Initialize first auxiliary vector
v1(2:2:end) = k_HB;            % Assign harmonic order to even indices
v2 = -v1;                      % Negative of v1 (for skew symmetry)

%% 2. Build harmonic operator matrices (Xi1: 1st-order, Xi2: 2nd-order)
% Xi1: Skew-symmetric matrix for first-order harmonic derivative
% Xi2: Square of Xi1 (second-order harmonic derivative)
Xi1 = diag(v1,1) + diag(v2,-1);
Xi2 = Xi1*Xi1;

%% 3. Generate sampling points for standard HB DFT (tau ∈ [0, 2π))
% taui: Uniformly spaced sampling points in time domain
taui = (0:2*pi/S_HB:2*pi-2*pi/S_HB)';

%% 4. Construct inverse DFT matrix (H_iDFT) and DFT matrix (H_DFT)
% H_iDFT: Inverse DFT matrix (time -> frequency), columns correspond to:
%         - Column 1: DC component (constant term)
%         - Even columns: Cosine terms (cos(k_HB*taui))
%         - Odd columns (≥3): Sine terms (sin(k_HB*taui))
H_iDFT = zeros(S_HB,L_HB);
H_iDFT(:,1) = ones(S_HB,1);                     % DC component (1s column)
H_iDFT(:,2:2:end) = cos(taui*k_HB);             % Cosine harmonic components
H_iDFT(:,3:2:end) = sin(taui*k_HB);             % Sine harmonic components

% H_DFT: DFT matrix (frequency -> time), scaled for orthogonality
%        - Transpose of H_iDFT (frequency domain dimension first)
%        - Scale non-DC terms by 2 (DFT orthogonality condition)
%        - Normalize by number of sampling points (S_HB)
H_DFT = H_iDFT';
H_DFT(2:end,:) = 2*H_DFT(2:end,:);              % Scale harmonic terms
H_DFT = H_DFT/(S_HB);                           % Normalization factor

%% 5. Generate weighted sampling points (higher-resolution for verification)
% W_taui: Uniformly spaced sampling points for weighted DFT (rad)
W_taui = (0:2*pi/SW_HB:2*pi-2*pi/SW_HB)';

%% 6. Construct weighted inverse DFT matrix (higher-resolution)
% W_H_iDFT: Weighted inverse DFT matrix (same structure as H_iDFT but more sampling points)
W_H_iDFT = zeros(SW_HB,L_HB);
W_H_iDFT(:,1) = ones(SW_HB,1);                  % DC component
W_H_iDFT(:,2:2:end) = cos(W_taui*k_HB);         % Cosine terms (weighted)
W_H_iDFT(:,3:2:end) = sin(W_taui*k_HB);         % Sine terms (weighted)

end

function [Xi1,Xi2,W_Gamma_0,W_Gamma_1] = FD(k_FD,U_FD,S_FD,SW_FD)
% FD - Finite Difference (FD) method for periodic derivative approximation
%      combined with periodic cubic spline interpolation
%
% INPUTS:
%   k_FD    - Array of grid point offsets (wavenumbers) for FD stencil
%             (defines which neighboring grid points are selected for FD derivative calculation,
%              e.g., k_FD = [-2,-1,0,1,2] means using 2 left, current, 2 right points)
%   U_FD    - Number of discrete points in one period [0, 2π)
%   S_FD    - Number of sampling points for standard FD calculation
%   SW_FD   - Number of high-resolution sampling points for spline interpolation
%
% OUTPUTS:
%   Xi1     - First-order periodic finite difference operator matrix (1/Δτ scaled)
%   Xi2     - Second-order periodic finite difference operator matrix (1/Δτ² scaled)
%   W_Gamma_0- Weighted interpolation matrix (periodic spline, SW_FD×U_FD)
%   W_Gamma_1- First-order derivative interpolation matrix (spline + Xi1)
%
% DESCRIPTION:
%   This function constructs periodic finite difference operators (1st/2nd order)
%   for harmonic components, then combines with periodic cubic spline interpolation
%   to generate high-resolution derivative matrices for periodic solution analysis.

%% 1. Input validation: U_FD and S_FD must have the same dimension
assert(numel(U_FD)==numel(S_FD), ...
    '%s:  U and S in FD is different');

%% 2. Calculate coefficients for FD operator (alp1: 1st order, alp2: 2nd order)
m = length(k_FD); % Dimension of harmonic wavenumber array (interval size)
% TY1/TY2: Target vectors for solving FD coefficient (1st/2nd derivative)
TY1 = zeros(m,1);
TY1(2) = 1;       % Target for first-order derivative coefficient
TY2 = zeros(m,1);
TY2(3) = 2*1;     % Target for second-order derivative coefficient

% K_FDM: Vandermonde matrix of wavenumbers (k_FD^0, k_FD^1, ..., k_FD^(m-1))
K_FDM = zeros(m,m);
for ii = 0:m-1
    K_FDM(ii+1,:) = k_FD.^ii;
end
% Solve linear system to get FD coefficients (alp1: 1st order, alp2: 2nd order)
alp1 = K_FDM\TY1;
alp2 = K_FDM\TY2;

% Numerical stabilization: set negligible coefficients to zero
for ii = 1:m
    if abs(alp1(ii))<10e-8
        alp1(ii) = 0;
    end
    if abs(alp2(ii))<10e-8
        alp2(ii) = 0;
    end
end

%% 3. Classify wavenumbers (zero/positive/negative) for periodic FD construction
[~,I_0] = find(k_FD==0);  % Index of zero wavenumber (DC component)
[~,I_p] = find(k_FD>0);  % Indices of positive wavenumbers
[~,I_n] = find(k_FD<0);  % Indices of negative wavenumbers

%% 4. Construct 1st/2nd order periodic FD operator matrices (Xi1/Xi2)
% Initialize with zero wavenumber component (diagonal matrix)
Xi1 = diag(ones(U_FD,1)*alp1(I_0),0);
Xi2 = diag(ones(U_FD,1)*alp2(I_0),0);

% Add positive wavenumber components (periodic off-diagonal terms)
for ii = 1:length(I_p)
    % First-order FD operator (positive wavenumber)
    Xi1 = Xi1 + diag(ones(U_FD-abs(k_FD(I_p(ii))),1)*alp1(I_p(ii)),k_FD(I_p(ii)));
    Xi1 = Xi1 + diag(ones(U_FD-(U_FD-abs(k_FD(I_p(ii)))),1)*alp1(I_p(ii)),-U_FD+k_FD(I_p(ii)));
    % Second-order FD operator (positive wavenumber)
    Xi2 = Xi2 + diag(ones(U_FD-abs(k_FD(I_p(ii))),1)*alp2(I_p(ii)),k_FD(I_p(ii)));
    Xi2 = Xi2 + diag(ones(U_FD-(U_FD-abs(k_FD(I_p(ii)))),1)*alp2(I_p(ii)),-U_FD+k_FD(I_p(ii)));
end

% Add negative wavenumber components (periodic off-diagonal terms)
for ii = 1:length(I_n)
    % First-order FD operator (negative wavenumber)
    Xi1 = Xi1 + diag(ones(U_FD-abs(k_FD(I_n(ii))),1)*alp1(I_n(ii)),k_FD(I_n(ii)));
    Xi1 = Xi1 + diag(ones(U_FD-(U_FD-abs(k_FD(I_n(ii)))),1)*alp1(I_n(ii)),U_FD+k_FD(I_n(ii)));
    % Second-order FD operator (negative wavenumber)
    Xi2 = Xi2 + diag(ones(U_FD-abs(k_FD(I_n(ii))),1)*alp2(I_n(ii)),k_FD(I_n(ii)));
    Xi2 = Xi2 + diag(ones(U_FD-(U_FD-abs(k_FD(I_n(ii)))),1)*alp2(I_n(ii)),U_FD+k_FD(I_n(ii)));
end

%% 5. Scale FD operators by time step (Δτ = 2π/U_FD)
Deltetau = 2*pi/U_FD; % Time step of periodic discrete points
Xi1 = Xi1/Deltetau;   % Scale 1st-order FD operator (d/dτ = 1/Δτ * FD)
Xi2 = Xi2/Deltetau^2; % Scale 2nd-order FD operator (d²/dτ² = 1/Δτ² * FD)

%% 6. Periodic cubic spline interpolation (high-resolution sampling)
[A,B,C,D] = periodic_spline(U_FD); % Get spline coefficient matrices
m = SW_FD/U_FD;                    % Interpolation refinement factor (SW_FD = m*U_FD)
SPLine = zeros(SW_FD,U_FD);        % Spline interpolation matrix (high-resolution)
h = Deltetau;                      % Original time step
SPLine(1:m:end,:) = eye(U_FD);     % Assign original discrete points (identity matrix)

% Construct spline interpolation matrix for refined points
for ii = 1:m-1
    deltexi = ii/m*h; % Normalized offset from original grid point
    % Cubic spline polynomial: S(τ) = A + B*Δτ + C*Δτ² + D*Δτ³
    SPLine(ii+1:m:end,:) = A + B*deltexi + C*deltexi^2 + D*deltexi^3;
end

%% 7. Generate weighted interpolation matrices (Gamma)
W_Gamma_0 = SPLine;        % Zero-order interpolation (high-resolution values)
W_Gamma_1 = SPLine*Xi1;    % First-order derivative interpolation (values × 1st FD)

end


function [A,B,C,D] = periodic_spline(U) % cubic spline interpolation for periodic solution
% PERIODIC_SPLINE - Construct periodic cubic spline interpolation matrices
%
% INPUT:
%   U    - Number of discrete points in one period [0, 2π)
%
% OUTPUTS:
%   A/B/C/D - Coefficient matrices for cubic spline polynomial:
%             S(x) = A·y + B·y' + C·y'' + D·y''' (y = discrete values)
%
% DESCRIPTION:
%   This subfunction builds periodic cubic spline matrices with cyclic boundary
%   conditions (continuous at τ=0 and τ=2π for all derivatives).

A = eye(U); % Identity matrix (zero-order term, y)
h = 2*pi/U; % Uniform grid spacing in [0, 2π)

%% Construct linear system for 2nd derivative (periodic boundary)
% C1: Coefficient matrix for 2nd derivative continuity (cyclic tridiagonal)
C1 = diag(ones(U,1)*4*h^2,0) + ...     % Main diagonal: 4h²
    diag(ones(U-1,1)*h^2,1) + ...     % Upper 1st diagonal: h²
    diag(ones(U-1,1)*h^2,-1) + ...    % Lower 1st diagonal: h²
    diag(h^2,U-1) + ...               % Cyclic term: (1,U) position
    diag(h^2,1-U);                    % Cyclic term: (U,1) position

% C2: Right-hand side matrix for 2nd derivative (cyclic boundary)
C2 = diag(-ones(U,1)*6,0) + ...        % Main diagonal: -6
    diag(ones(U-1,1)*3,1) + ...       % Upper 1st diagonal: 3
    diag(ones(U-1,1)*3,-1) + ...      % Lower 1st diagonal: 3
    diag(3,U-1) + ...                 % Cyclic term: (1,U) position
    diag(3,1-U);                      % Cyclic term: (U,1) position

C =  C1\C2; % Solve linear system C1·C = C2 (2nd derivative coefficients)

%% Construct 3rd derivative coefficient matrix (D)
D1 = diag(-ones(U,1)/h/3,0) + ...      % Main diagonal: -1/(3h)
    diag(ones(U-1,1)/h/3,1) + ...     % Upper 1st diagonal: 1/(3h)
    diag(1/h/3,1-U);                  % Cyclic term: (1,U) position
D = D1*C; % 3rd derivative coefficients (from 2nd derivative)

%% Construct 1st derivative coefficient matrix (B)
B1 = diag(-ones(U,1)/h,0) + ...        % Main diagonal: -1/h
    diag(ones(U-1,1)/h,1) + ...       % Upper 1st diagonal: 1/h
    diag(1/h,1-U);                    % Cyclic term: (1,U) position
% Combine terms for 1st derivative (B = B1 - hC - h²D)
B = B1 + (-h)*C + (-h^2)*D;
end


function [Xi_i,W_Xi_i,W_Gamma_0,W_Gamma_1] = CO(k_CO,U_CO,S_CO,SW_CO)
% CO - Collocation Method (CO) for periodic solution approximation
%      Supports Gauss collocation points, base collocation points, and high-resolution extension
%
% INPUTS:
%   k_CO    - 2-element array defining collocation configuration:
%             k_CO(1) = Degree of collocation polynomial (m)
%             k_CO(2) = Collocation point type flag:
%                       1 = Gauss collocation points (high precision)
%                       2 = Base collocation points (uniform grid points)
%   U_CO    - Total number of discrete base points in one period [0, 2π)
%   S_CO    - Number of sampling points for standard collocation calculation (match U_CO)
%   SW_CO   - Number of high-resolution sampling points for collocation extension
%
% OUTPUTS:
%   Xi_i    - Cell array of collocation point operators (3×1) 
%             Xi_i{1} = 0th-order collocation point matrix (value interpolation at collocation points)
%             Xi_i{2} = 1st-order collocation point operator (1st derivative at collocation points)
%             Xi_i{3} = 2nd-order collocation point operator (2nd derivative at collocation points)
%   W_Xi_i  - Cell array of base point operators (2×1) 
%             W_Xi_i{1} = 1st-order base point operator (1st derivative at base points)
%             W_Xi_i{2} = 2nd-order base point operator (2nd derivative at base points)
%   W_Gamma_0- High-resolution 0th-order interpolation matrix (periodic collocation extension)
%   W_Gamma_1- High-resolution 1st-order derivative interpolation matrix (periodic collocation extension)
%
% DESCRIPTION:
%   This function constructs collocation point operators and base point operators for periodic solution analysis,
%   including Gauss/base collocation point schemes and high-resolution collocation extension.
%   Periodic boundary conditions are enforced via cyclic matrix construction for all operators.

%% 1. Input validation: Ensure U_CO and S_CO have matching dimensions
assert(numel(U_CO)==numel(S_CO), ...
    '%s:  U and S in CO is different');

%% 2. Initialize output cell arrays
Xi_i = cell(3,1);   % Collocation point operators (0th/1st/2nd order) 
W_Xi_i = cell(2,1); % Base point operators (1st/2nd order) 

%% 3. Parse collocation configuration parameters
m = k_CO(1);     % Degree of the collocation polynomial (points per collocation interval)
flag = k_CO(2);  % Collocation point type flag: 1=Gauss points, 2=base points
Deltetau = 2*pi/U_CO;  % Uniform time step of periodic base points [0, 2π)
N = U_CO/m;            % Number of collocation intervals (total base points / polynomial degree)

%% 4. Collocation & base point operator construction (branch by point type)
if flag == 1
    % --------------------------
    % Case 1: Gauss collocation points (high-precision, collocation≠base)
    % --------------------------
    % Generate and normalize Gauss collocation nodes to [0, m] interval
    c_node = gauss_node(m);                  % Raw Gauss nodes ([-1, 1])
    c_node = m/2*(c_node+1);                % Rescale to [0, m] (match collocation interval)
    base_p = (0:m-1)';                       % Uniform base points indices [0,1,...,m-1]
    
    % Construct matrices for Lagrange polynomial denominator calculation
    a1 = repmat( (0:m)',1,m+1 );             % Repeat polynomial degrees
    a1 = a1(:);                              % Vectorize matrix
    n = m+1;
    a1(1:n+1:n*n) = [];                      % Remove diagonal elements (avoid division by zero)
    a1 = reshape(a1,n-1,n);                  % Reshape to (m)×(m+1) matrix
    a2 = repmat( (0:m),m,1 );                % Repeat base points for denominator
    
    % Calculate denominator/numerator of Lagrange basis polynomials
    lo_denominator = prod( a2-a1 );          % Denominator of Lagrange polynomial
    molecule_w = kron(c_node',ones(m,m+1)) - kron( ones(1,m),a1 );  % Numerator (collocation points)
    molecule_b = kron(base_p',ones(m,m+1)) - kron( ones(1,m),a1 );  % Numerator (base points)
    
    % 0th-order Lagrange basis (value interpolation at collocation points)
    l0_molecule_w = prod(molecule_w);                            % Product of numerator terms
    l0_w = reshape( l0_molecule_w./repmat( lo_denominator,1,m ), m+1, m)';  % Normalize and reshape
    
    % 1st-order Lagrange basis (1st derivative for collocation/base points)
    l1_molecule_w = zeros(size(l0_molecule_w));  % Collocation point derivative numerator
    l1_molecule_b = zeros(size(l0_molecule_w));  % Base point derivative numerator
    for ii = 1:m
        % Remove ii-th row for derivative calculation (Lagrange derivative formula)
        a4 = molecule_w;
        a4(ii,:) = [];
        l1_molecule_w = prod(a4) + l1_molecule_w;  % Collocation point 1st derivative
        
        a5 = molecule_b;
        a5(ii,:) = [];
        l1_molecule_b = prod(a5) + l1_molecule_b;  % Base point 1st derivative
    end
    % Reshape and scale by time step (d/dτ = 1/Δτ * discrete derivative)
    l1_w = reshape( l1_molecule_w./repmat( lo_denominator,1,m ), m+1, m)'/Deltetau; % Collocation point 1st-order
    l1_b = reshape( l1_molecule_b./repmat( lo_denominator,1,m ), m+1, m)'/Deltetau; % Base point 1st-order
    
    % Construct periodic 0th-order collocation point matrix
    Lag0_w_a = kron(eye(N),l0_w(:,1:m));     % Block diagonal matrix (main interval)
    Lag0_w_b = kron(eye(N),[zeros(m,m-1),l0_w(:,1+m)]);  % Periodic boundary block
    Lag0_w_b = [Lag0_w_b(:,end),Lag0_w_b(:,1:end-1)];     % Cyclic shift for periodicity
    Lag0_w = Lag0_w_a + Lag0_w_b;            % 0th-order operator at collocation points
    
    % Construct periodic 1st-order collocation point operator
    Lag1_w_a = kron(eye(N),l1_w(:,1:m));     % Block diagonal matrix (main interval)
    Lag1_w_b = kron(eye(N),[zeros(m,m-1),l1_w(:,1+m)]);  % Periodic boundary block
    Lag1_w_b = [Lag1_w_b(:,end),Lag1_w_b(:,1:end-1)];     % Cyclic shift for periodicity
    Lag1_w = Lag1_w_a + Lag1_w_b;            % 1st-order operator at collocation points
    
    % Construct periodic 1st-order base point operator
    Lag1_b_a = kron(eye(N),l1_b(:,1:m));     % Block diagonal matrix (main interval)
    Lag1_b_b = kron(eye(N),[zeros(m,m-1),l1_b(:,1+m)]);  % Periodic boundary block
    Lag1_b_b = [Lag1_b_b(:,end),Lag1_b_b(:,1:end-1)];     % Cyclic shift for periodicity
    Lag1_b = Lag1_b_a + Lag1_b_b;            % 1st-order operator at base points

else
    % --------------------------
    % Case 2: Base collocation points (uniform grid, collocation=base)
    % --------------------------
    c_node = (0:m-1)';                       % Collocation points = uniform base points
    % Construct matrices for Lagrange polynomial denominator calculation
    a1 = repmat( (0:m)',1,m+1 );             % Repeat polynomial degrees
    a1 = a1(:);                              % Vectorize matrix
    n = m+1;
    a1(1:n+1:n*n) = [];                      % Remove diagonal elements (avoid division by zero)
    a1 = reshape(a1,n-1,n);                  % Reshape to (m)×(m+1) matrix
    a2 = repmat( (0:m),m,1 );                % Repeat base points for denominator
    
    % Calculate denominator/numerator of Lagrange basis polynomials
    lo_denominator = prod( a2-a1 );          % Denominator of Lagrange polynomial
    molecule_w = kron(c_node',ones(m,m+1)) - kron( ones(1,m),a1 );  % Numerator (collocation=base)
    
    % 1st-order Lagrange basis (1st derivative at collocation/base points)
    l1_molecule_w = zeros(1,size(molecule_w,2));
    for ii = 1:m
        a4 = molecule_w;
        a4(ii,:) = [];
        l1_molecule_w = prod(a4) + l1_molecule_w;  % 1st derivative (collocation=base)
    end
    % Reshape and scale by time step
    l1_w = reshape( l1_molecule_w./repmat( lo_denominator,1,m ), m+1, m)'/Deltetau; % 1st-order operator (collocation=base)
    
    % 0th-order collocation point matrix (identity, collocation=base)
    Lag0_w = eye(U_CO);                      % Exact value at collocation/base points
    
    % Construct periodic 1st-order collocation point operator (collocation=base)
    Lag1_w_a = kron(eye(N),l1_w(:,1:m));     % Block diagonal matrix (main interval)
    Lag1_w_b = kron(eye(N),[zeros(m,m-1),l1_w(:,1+m)]);  % Periodic boundary block
    Lag1_w_b = [Lag1_w_b(:,end),Lag1_w_b(:,1:end-1)];     % Cyclic shift for periodicity
    Lag1_w = Lag1_w_a + Lag1_w_b;            % 1st-order operator at collocation points
    
    % 1st-order base point operator (same as collocation, collocation=base)
    Lag1_b = Lag1_w;
end

%% 5. Assign operators to output cell arrays
% Xi_i: Collocation point operators 
Xi_i{1} = Lag0_w;                  % 0th-order (value interpolation at collocation points)
Xi_i{2} = Lag1_w;                  % 1st-order (1st derivative at collocation points)
Xi_i{3} = Lag1_w*Lag1_b;           % 2nd-order (2nd derivative at collocation points)
% W_Xi_i: Base point operators 
W_Xi_i{1} = Lag1_b;                % 1st-order (1st derivative at base points)
W_Xi_i{2} = Lag1_b*Lag1_b;         % 2nd-order (2nd derivative at base points)

%% 6. High-resolution collocation extension (periodic interpolation)
% Calculate high-resolution collocation parameters
mm_hr = SW_CO/N;                     % High-resolution points per collocation interval
c_node_hr = (0:mm_hr-1)'/mm_hr*m;    % High-resolution nodes (normalized to [0,m])

% Reconstruct Lagrange polynomial denominator (reuse base logic)
a1 = repmat( (0:m)',1,m+1 );
a1 = a1(:);
n = m+1;
a1(1:n+1:n*n) = [];
a1 = reshape(a1,n-1,n);
a2 = repmat( (0:m),m,1 );
lo_denominator = prod( a2-a1 );

% Calculate numerator for high-resolution Lagrange basis
molecule_w_hr = kron(c_node_hr',ones(m,m+1)) - kron( ones(1,mm_hr),a1 );  % High-resolution numerator

% 0th-order high-resolution Lagrange basis
l0_molecule_w_hr = prod(molecule_w_hr);                                  % Product of numerator terms
l0_w_hr = reshape( l0_molecule_w_hr./repmat( lo_denominator,1,mm_hr ), m+1, mm_hr)';  % Normalize/reshape

% 1st-order high-resolution Lagrange basis (1st derivative)
l1_molecule_w_hr = zeros(1,size(molecule_w_hr,2));
for ii = 1:m
    a4 = molecule_w_hr;
    a4(ii,:) = [];
    l1_molecule_w_hr = prod(a4) + l1_molecule_w_hr;  % High-resolution 1st derivative
end
% Reshape and scale high-resolution derivative matrix
l1_w_hr = reshape( l1_molecule_w_hr./repmat( lo_denominator,1,mm_hr ), m+1, mm_hr)'/Deltetau;

% Construct periodic 0th-order high-resolution interpolation matrix
Lag0_w_a_hr = kron(eye(N),l0_w_hr(:,1:m));                     % Block diagonal matrix (main interval)
Lag0_w_b_hr = kron(eye(N),[zeros(mm_hr,m-1),l0_w_hr(:,1+m)]);  % Periodic boundary block
Lag0_w_b_hr = [Lag0_w_b_hr(:,end),Lag0_w_b_hr(:,1:end-1)];     % Cyclic shift for periodicity
W_Gamma_0 = Lag0_w_a_hr + Lag0_w_b_hr;                         % High-resolution 0th-order interpolation

% Construct periodic 1st-order high-resolution derivative matrix
Lag1_w_a_hr = kron(eye(N),l1_w_hr(:,1:m));                     % Block diagonal matrix (main interval)
Lag1_w_b_hr = kron(eye(N),[zeros(mm_hr,m-1),l1_w_hr(:,1+m)]);  % Periodic boundary block
Lag1_w_b_hr = [Lag1_w_b_hr(:,end),Lag1_w_b_hr(:,1:end-1)];     % Cyclic shift for periodicity
W_Gamma_1 = Lag1_w_a_hr + Lag1_w_b_hr;                         % High-resolution 1st-order derivative

end

%% Subfunction: Generate Gauss collocation nodes (symbolic calculation)
function c_node = gauss_node(m)
% GAUSS_NODE - Generate Gauss collocation nodes via symbolic differentiation
%
% INPUT:
%   m    - Degree of the collocation polynomial
%
% OUTPUT:
%   c_node - Sorted Gauss collocation nodes in [-1, 1] (roots of m-th Legendre polynomial)
%
% DESCRIPTION:
%   Calculates Gauss nodes by solving the m-th derivative of the Legendre polynomial ((t^2-1)^m),
%   which are the optimal collocation points for high-precision numerical integration/derivation.

g_l = ((sym('t'))^2-1)^m;       % Legendre polynomial (symbolic expression)
dg_l = diff(g_l,m);             % m-th derivative (Gauss node condition: roots of m-th Legendre polynomial)
c_node = double(solve(dg_l));   % Solve for symbolic roots -> numerical Gauss nodes
c_node = sort(c_node);          % Sort nodes in ascending order [-1, 1]
end
