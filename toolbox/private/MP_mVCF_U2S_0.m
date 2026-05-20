function X0 = MP_mVCF_U2S_0(PHI,Xd,U_prod)
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
Xii = Xd;
% Iterative spectral projection (backward over extended frequencies)
for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    % Reshape + transpose for spectral transformation
    Xii_wave = pagetranspose( reshape(Xii,pSiip1Td,U(kk),U1Tii) );
    Xii_wave = reshape(Xii_wave,U(kk),pSiip1Td*U1Tii);
    % Spectral method (HB/CO/FD)
    switch PHI.method{kk}
        case{'HB'} % Harmonic Balance (HB)
            Xiiminus1_bar = PHI.Gam0{kk}*Xii_wave;
        case{'CO'} % Collocation (CO) 
            Xiiminus1_bar = Xii_wave;                     
        case{'FD'} % Finite Difference (FD)
            Xiiminus1_bar = Xii_wave;
    end
    pSiip1Td = pSiip1Td*pS(kk);
    Xii = reshape( Xiiminus1_bar, pSiip1Td, U1Tii  );
end
X0 = Xii';
end