function [X0,dX0] = MP_mVCF_U2S_01(PHI,Xd,Om,U_prod)
%MP_MVCF_U2S_01  Reconstruct state and first hyper-time derivative samples.
%
% The routine evaluates the state represented by Xd and its derivative with
% respect to each hyper-time direction. The derivative matrices stored in PHI
% are combined with the current base frequencies Om.
%
% On input:
%   PHI    : m-VCF basis and derivative matrices.
%   Xd     : U-domain coefficient row for one state component.
%   Om     : Current base-frequency vector.
%   U_prod : Total coefficient count.
%
% On output:
%   X0     : State samples.
%   dX0    : Frequency-scaled derivative samples used by F.
pSiip1Td = 1;
U1Tii = U_prod;
U = cell2mat(PHI.U_omod);
pS = cell2mat(PHI.pS_omod);
Xii  = Xd;
dXii = zeros(size(Xd));
% Iterative spectral projection (backward over extended frequencies)
for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    % Reshape + transpose for spectral transformation
    Xii_wave  = pagetranspose( reshape(Xii,pSiip1Td,U(kk),U1Tii) );
    Xii_wave  = reshape(Xii_wave,U(kk),pSiip1Td*U1Tii);
    dXii_wave = pagetranspose( reshape(dXii,pSiip1Td,U(kk),U1Tii) );
    dXii_wave = reshape(dXii_wave,U(kk),pSiip1Td*U1Tii);
    % Spectral method (HB/CO/FD)
    switch PHI.me_omod{kk}
        case{'HB'} % Harmonic Balance (HB)
            Xiiminus1_bar  = PHI.Gam0{kk}*Xii_wave;
            dXiiminus1_bar = Om(kk)*PHI.Gam1{kk}*Xii_wave+...
                PHI.Gam0{kk}*dXii_wave;
        case{'CO'} % Collocation (CO) 
            Xiiminus1_bar  = Xii_wave;
            dXiiminus1_bar = Om(kk)*PHI.Gam1{kk}*Xii_wave+...
                dXii_wave;
        case{'FD'} % Finite Difference (FD)
            Xiiminus1_bar  = Xii_wave;
            dXiiminus1_bar = Om(kk)*PHI.Gam1{kk}*Xii_wave+...
                dXii_wave;
    end
    pSiip1Td = pSiip1Td*pS(kk);
    Xii  = reshape( Xiiminus1_bar, pSiip1Td, U1Tii  );
    dXii = reshape( dXiiminus1_bar, pSiip1Td, U1Tii  );    
end
X0  = Xii';
dX0 = dXii';
end