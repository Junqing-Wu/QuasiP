function Nd = MP_mVCF_S2U(PHI,N0,pS_prod)
%MP_MVCF_S2U  Transform nonlinear samples from the S-domain to U-domain.
%
% Nonlinear forces are evaluated at PHI.pS_omod sample points. This routine
% applies the inverse basis transform direction by direction, so the sampled
% values N0 are represented by the coefficient vector Nd used in the m-VCF
% residual. The order of the tensor-product coefficients follows PHI.U_omod.
%
% On input:
%   PHI    : m-VCF discretization data with U, pS, and iGam0 transforms.
%   N0     : Sample-domain values arranged along the product sample grid.
%   pS_prod: Product of the sample counts in all hyper-time directions.
%
% On output:
%   Nd     : U-domain coefficient vector in the ordering used by Xd.
% Tensor-product projection used by the code:
%   Nd = S2U(PHI,N0), where N0 contains samples on the S-domain grid.
% Nd = S2U(PHI,N0) maps nonlinear samples N0 on the S-domain grid to U-domain coefficients.
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