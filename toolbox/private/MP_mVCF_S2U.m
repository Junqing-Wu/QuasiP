function Nd = MP_mVCF_S2U(PHI,N0,pS_prod)
% MP_mVCF_S2U: S-domain (spectral) → U-domain (harmonic) mapping for
%              Nonlinear terms or Excitation
% Inputs:
%   - PHI: MP parametrization config (HB/CO/FD methods, U/S dimensions)
%   - N0: S-domain nonlinear force vector (dimension: 1 × pS_prod )
%   - S_prod: Total S-domain spectral combinations
% Output:
%   - Nd: U-domain nonlinear force matrix (dimension: U_prod × 1)
U = cell2mat(PHI.U_omod);
pS = cell2mat(PHI.pS_omod);
pS_iiTd = pS_prod;
U_1Tiiminus1 = 1;
Niiminus1 = N0;
% Iterative inverse spectral projection (forward over extended frequencies)
for kk = 1:PHI.d
    pS_iiTd = pS_iiTd/pS(kk);
    % Reshape + transpose for inverse spectral transformation
    Niiminus1_bar = pagetranspose(reshape(Niiminus1,U_1Tiiminus1,pS(kk),pS_iiTd) );
    Niiminus1_bar = reshape(Niiminus1_bar,pS(kk),U_1Tiiminus1*pS_iiTd);

    % Inverse spectral method (HB/CO/FD)
    switch PHI.me_omod{kk}
        case{'HB'} % Inverse Harmonic Balance
            Fnlii = PHI.iGam0{kk}*Niiminus1_bar;
        case{'CO'} % Collocation (identity mapping)
            Fnlii = Niiminus1_bar;
        case{'FD'} % Finite Difference (identity mapping)
            Fnlii = Niiminus1_bar;
    end
    U_1Tiiminus1 = U_1Tiiminus1*U(kk);
    Niiminus1 = reshape(Fnlii,U_1Tiiminus1,pS_iiTd);
end
Nd = Niiminus1;
end
