function [X0,X0dot] = MP_mVCF_U2S_xv(Xd,Om,PHI,U_prod)
pwSiip1Td = 1;
U1Tii = U_prod;
U = cell2mat(PHI.U_omod);
pwS = cell2mat(PHI.pwS_omod);

Xii = reshape(Xd,pwSiip1Td,U1Tii);
Xiidot = zeros(size(Xii));
for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    Xii_wave = pagetranspose( reshape(Xii,pwSiip1Td,U(kk),U1Tii) );
    Xii_wave = reshape(Xii_wave,U(kk),pwSiip1Td*U1Tii);

    Xiidot_wave = pagetranspose( reshape(Xiidot,pwSiip1Td,U(kk),U1Tii) );
    Xiidot_wave = reshape(Xiidot_wave,U(kk),pwSiip1Td*U1Tii);
   
    GAM0 = PHI.wGam0{kk};
    GAM1 = PHI.wGam1{kk};

    Xiiminus1_bar = GAM0*Xii_wave;
    Xiidotminus1_bar = Om(kk)*GAM1*Xii_wave+...
                       GAM0*Xiidot_wave;

    pwSiip1Td = pwSiip1Td*pwS(kk);
    Xii = reshape( Xiiminus1_bar, pwSiip1Td, U1Tii  );
    Xiidot = reshape( Xiidotminus1_bar, pwSiip1Td, U1Tii  );
end
X0 = Xii;
X0dot = Xiidot;
end


