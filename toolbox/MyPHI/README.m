%% MyPHI is used to solve for the invariance matrix involved in the construction of numerical methods. 

%% 1: parametrization_tori: MP (Multi-periods parametrization for tori)

% 1.1: MP.type: mVCF (multi-steps variable-coefficients formulation)

% 1.1.1 MP.mVCF.method: HB(harmonic balance method)
% 1.1.1.1 MP.mVCF.HB  ---> k:(1:1:h)                   harmonic order  
%                     ---> U: 2*numel(k)+1             number of Phi_i functions
%                     ---> S:                          number of tau_i_bar for computation of tori
%                     ---> S_wave:                     number of tau_i_bar for post-processing of tori

% 1.1.2 MP.mVCF.method: CO(collocation method)
% 1.1.2.1 MP.mVCF.HB  ---> k:[m,isgauss]               m: order of Lagrange polynomial function 
%                     ---> U:                          number of Phi_i functions (U = mP)
%                     ---> S:(U=S)                     number of tau_i_bar for computation of tori
%                     ---> S_wave:(S_w/S = 1,2...)     number of tau_i_bar for post-processing of tori

% 1.1.3 MP.mVCF.method: FD(finite difference method)
% 1.1.3.1 MP.mVCF.HB  ---> k:-w_l:1:w_r                W = {-w_l,...,w_r} 
%                     ---> U:                          number of Phi_i functions
%                     ---> S:(U=S)                     number of tau_i_bar for computation of tori
%                     ---> S_wave:(S_w/S = 1,2...)     number of tau_i_bar for post-processing of tori

% 1.2: MP.type: sCCF (single-step constant-coefficients formulation)


%% 2: parametrization_tori: MT (Multi-trajectories parametrization for tori)


% 2.1: MT.type: DM (discretization method)


% 2.1: MT.type: SHOOT (shooting method)