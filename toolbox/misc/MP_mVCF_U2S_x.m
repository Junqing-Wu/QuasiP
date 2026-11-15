function X0 = MP_mVCF_U2S_x(Xd,PHI,U_prod)
%MP_MVCF_U2S_X  Reconstruct m-VCF state samples for Processor routines.
%
% Apply the tensor-product PHI.wGam0 transforms direction by direction to
% evaluate the physical state on the dense S_wave grid.
%
% On input:
%   Xd      : Coefficients of one state component in U-domain ordering.
%   PHI     : m-VCF basis, transform, period, and sample-count data.
%   U_prod  : Product of the U_i coefficient counts.
%
% On output:
%   X0      : State samples on the dense hyper-time grid.

pwSiip1Td = 1;
U1Tii = U_prod;
U = cell2mat(PHI.U_omod);
pwS = cell2mat(PHI.pwS_omod);

Xii = reshape(Xd,pwSiip1Td,U1Tii);

for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    Xii_wave = pagetranspose( reshape(Xii,pwSiip1Td,U(kk),U1Tii) );
    Xii_wave = reshape(Xii_wave,U(kk),pwSiip1Td*U1Tii);
   
    GAM0 = PHI.wGam0{kk};
    Xiiminus1_bar = GAM0*Xii_wave;

    pwSiip1Td = pwSiip1Td*pwS(kk);
    Xii = reshape( Xiiminus1_bar, pwSiip1Td, U1Tii  );
end

X0 = Xii;
end