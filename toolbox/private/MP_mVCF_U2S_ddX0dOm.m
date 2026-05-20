function ddX0dOm = MP_mVCF_U2S_ddX0dOm(PHI,Xd,U_prod,dd)
% MP_mVCF_U2S_0: U-domain (harmonic) → S-domain (spectral) mapping for QP
% Inputs:
%   - PHI: MP parametrization config (HB/CO/FD methods, U/S dimensions)
%   - Xd: U-domain coefficient vector (dimension: U_prod)
%   - U_prod: Total U-domain harmonic combinations
% Output:
%   - X0: S-domain spectral vector (dimension: S_prod)
pSiip1Td = 1;
U1Tii = U_prod;
U = cell2mat(PHI.U);
pS = cell2mat(PHI.pS);
ddXiidOm = Xd;
% Iterative spectral projection (backward over extended frequencies)
for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    % Reshape + transpose for spectral transformation
    ddXiidOm_wave = pagetranspose( reshape(ddXiidOm,pSiip1Td,U(kk),U1Tii) );
    ddXiidOm_wave = reshape(ddXiidOm_wave,U(kk),pSiip1Td*U1Tii);
    % Spectral method (HB/CO/FD)
    switch PHI.method{kk}
        case{'HB'} % Harmonic Balance (HB)
            if dd == kk
                ddXiiminus1_bardOm = PHI.Gam1{kk}*ddXiidOm_wave;
            else
                ddXiiminus1_bardOm = PHI.Gam0{kk}*ddXiidOm_wave;
            end
        case{'CO'} % Collocation (CO)
            if dd == kk
                ddXiiminus1_bardOm = PHI.Gam1{kk}*ddXiidOm_wave;
            else
                ddXiiminus1_bardOm = ddXiidOm_wave;
            end
        case{'FD'} % Finite Difference (FD)
            if dd == kk
                ddXiiminus1_bardOm = PHI.Gam1{kk}*ddXiidOm_wave;
            else
                ddXiiminus1_bardOm = ddXiidOm_wave;
            end
    end
    pSiip1Td = pSiip1Td*pS(kk);
    ddXiidOm = reshape( ddXiiminus1_bardOm, pSiip1Td, U1Tii  );
end
ddX0dOm = ddXiidOm;
end