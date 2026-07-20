function Nd_Xd = MP_mVCF_S2U_withXDX(PHI, N0_X0, N0_dX0, Om, pS_prod)
% MP_mVCF_S2U_withXDX: S-domain → U-domain mapping for nonlinear force gradients (∂N/∂Xd)
% Inputs:
%   - PHI: MP parametrization config (HB/CO/FD methods, U/S dimensions)
%   - N0_X0: S-domain gradient tensor (dimension: 1×pS_prod)
%   - S_prod: Total S-domain spectral combinations
% Output:
%   - Nd_Xd: U-domain gradient matrix (dimension: U_prod × U_prod)
U            = cell2mat(PHI.U_omod);
pS           = cell2mat(PHI.pS_omod);
pS_iiTd      = pS_prod;
U_1Tiiminus1 = 1;
N_Xiiminus1   = reshape(N0_X0,1,[]);
N_dXiiminus1  = reshape(N0_dX0,1,[]);
% Iterative gradient projection (forward over extended frequencies)
for kk = 1:PHI.d
    pS_iiTd = pS_iiTd/pS(kk);
    % Reshape + transpose for gradient transformation
    N_Xiiminus1_bar  = pagetranspose( reshape(N_Xiiminus1,U_1Tiiminus1^2,pS(kk),pS_iiTd) );
    N_Xiiminus1_bar  = reshape( N_Xiiminus1_bar, pS(kk),1,U_1Tiiminus1^2*pS_iiTd );
    N_dXiiminus1_bar = pagetranspose( reshape(N_dXiiminus1,U_1Tiiminus1^2,pS(kk),pS_iiTd) );
    N_dXiiminus1_bar = reshape( N_dXiiminus1_bar, pS(kk),1,U_1Tiiminus1^2*pS_iiTd );
    % Gradient transformation (HB/CO/FD)
    switch PHI.me_omod{kk}
        case{'HB'} % Harmonic Balance gradient
            N_Xiiminus1_bar_SU = repmat(N_Xiiminus1_bar,1,U(kk));
            N_dXiiminus1_bar_SU = repmat(N_dXiiminus1_bar,1,U(kk));

            N_Xii  = N_Xiiminus1_bar_SU.*PHI.Gam0{kk} +...
                     Om(kk)*N_dXiiminus1_bar_SU.*PHI.Gam1{kk};
            N_Xii  = pagemtimes( PHI.iGam0{kk},N_Xii );

            N_dXii = N_dXiiminus1_bar_SU.*PHI.Gam0{kk};
            N_dXii = pagemtimes( PHI.iGam0{kk},N_dXii );
        case{'CO'} % Collocation gradient (diagonal mapping)
            N_Xiiminus1_bar_SS  = zeros(pS(kk),pS(kk),U_1Tiiminus1^2*pS_iiTd);
            N_dXiiminus1_bar_SU = repmat(N_dXiiminus1_bar,1,U(kk));
            N_dXiiminus1_bar_SS = zeros(pS(kk),pS(kk),U_1Tiiminus1^2*pS_iiTd);
          
            for jj = 1:U_1Tiiminus1^2*pS_iiTd
                N_Xiiminus1_bar_SS(:,:,jj)  = diag(N_Xiiminus1_bar(:,:,jj));
                N_dXiiminus1_bar_SS(:,:,jj) = diag(N_dXiiminus1_bar(:,:,jj));                
            end

            N_Xii = N_Xiiminus1_bar_SS +...
                    Om(kk)*N_dXiiminus1_bar_SU.*PHI.Gam1{kk};
            N_dXii = N_dXiiminus1_bar_SS;
        case{'FD'} % Finite Difference gradient (diagonal mapping)
            N_Xiiminus1_bar_SS  = zeros(pS(kk),pS(kk),U_1Tiiminus1^2*pS_iiTd);
            N_dXiiminus1_bar_SU = repmat(N_dXiiminus1_bar,1,U(kk));
            N_dXiiminus1_bar_SS = zeros(pS(kk),pS(kk),U_1Tiiminus1^2*pS_iiTd);
          
            for jj = 1:U_1Tiiminus1^2*pS_iiTd
                N_Xiiminus1_bar_SS(:,:,jj)  = diag(N_Xiiminus1_bar(:,:,jj));
                N_dXiiminus1_bar_SS(:,:,jj) = diag(N_dXiiminus1_bar(:,:,jj));                
            end

            N_Xii = N_Xiiminus1_bar_SS +...
                    Om(kk)*N_dXiiminus1_bar_SU.*PHI.Gam1{kk};
            N_dXii = N_dXiiminus1_bar_SS;
    end
    % Reshape + transpose for next iteration
    N_Xiiminus1  = pagetranspose( reshape(N_Xii,U(kk),U(kk)*U_1Tiiminus1,U_1Tiiminus1,pS_iiTd) );
    N_dXiiminus1 = pagetranspose( reshape(N_dXii,U(kk),U(kk)*U_1Tiiminus1,U_1Tiiminus1,pS_iiTd) );
    
    U_1Tiiminus1 = U_1Tiiminus1*U(kk);

    N_Xiiminus1  = reshape(N_Xiiminus1,U_1Tiiminus1^2,pS_iiTd);
    N_dXiiminus1 = reshape(N_dXiiminus1,U_1Tiiminus1^2,pS_iiTd);    
end
Nd_Xd = transpose( reshape(N_Xiiminus1,U_1Tiiminus1,U_1Tiiminus1) );
end

