function X0 = MP_mVCF_U2S_0(PHI,Xd,U_prod)
%MP_MVCF_U2S_0  Reconstruct physical samples from m-VCF coefficients.
%
% Apply the tensor-product basis transforms PHI.Gam0 to convert one state
% component from U-domain coefficients to sample-domain values. This is used
% before evaluating the nonlinear force F at the hyper-time sample points.
%
% On input:
%   PHI    : m-VCF basis and transform data.
%   Xd     : U-domain coefficient row for one state component.
%   U_prod : Total number of U-domain coefficients.
%
% On output:
%   X0     : Reconstructed state samples.
% Tensor-product reconstruction used by the code:
%   X0 = U2S(PHI,Xd), with X0 sampled on the S-domain grid.
% X0 = U2S(PHI,Xd) reconstructs the physical samples from coefficient vector Xd.
pSiip1Td = 1;
U1Tii = U_prod;
U = cell2mat(PHI.U_omod);
pS = cell2mat(PHI.pS_omod);
Xii = Xd;
% Iterative spectral projection (backward over extended frequencies)
for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    % Reshape + transpose for spectral transformation
    Xii_wave = pagetranspose( reshape(Xii,pSiip1Td,U(kk),U1Tii) );
    Xii_wave = reshape(Xii_wave,U(kk),pSiip1Td*U1Tii);
    % Spectral method (HB/CO/FD)
    switch PHI.me_omod{kk}
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