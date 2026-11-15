function ddX0dOm = MP_mVCF_U2S_ddX0dOm(PHI,Xd,U_prod,dd)
% Frequency sensitivity:
%   dX0/dOm(dd) = d[U2S(PHI,Xd)]/dOm(dd).
% Only the derivative operator associated with direction dd contributes.
%MP_MVCF_U2S_DDX0DOM  Compute the derivative of reconstructed samples with respect to Om.
%
% Differentiate the U-to-S reconstruction in the requested frequency
% direction dd. HB uses PHI.Gam1; CO and FD use their corresponding first-
% derivative interpolation operator. The result is used when differentiating
% the m-VCF residual with respect to base frequencies.
%
% On input:
%   PHI : m-VCF basis and derivative data.
%   Xd : m-VCF coefficient vector.
%   U_prod : Product of coefficient counts.
%   dd : Frequency-direction index.
%                      frequency-direction index.
%
% On output:
%   ddX0dOm : Derivative of reconstructed state samples with respect to Om(dd).
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