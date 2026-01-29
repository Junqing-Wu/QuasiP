function C = MP_mVCF_INITIAL_C_byPHI(data,X_MF,OM)

PHI = data.PHI;
n = data.DOF;
S = PHI.S;
U = PHI.U;
d = PHI.d;
tau = PHI.tau;

k_h = X_MF(:,1:d)';
Q   = X_MF(:,1+d:end);

H_IDF = exp(1i*tau*k_h);
Xiiminus1 = real(H_IDF*Q);
% compute discrete Zii
S_iiTd = prod(cell2mat(S));
U_1Tiiminus1 = 1;
Xiiminus1 = reshape(Xiiminus1,U_1Tiiminus1,S_iiTd,n);
for ii = 1:d
    S_iiTd = S_iiTd/S{ii};
    Xiiminus1_bar = pagetranspose(reshape(Xiiminus1,U_1Tiiminus1,S{ii},S_iiTd,n) );
    iGAM = PHI.Gamma_i{3,ii};
    switch PHI.method{ii}
        case{'HB'}
            Zii = pagemtimes(iGAM,Xiiminus1_bar);
        case{'CO'}
            Zii = Xiiminus1_bar;
        case{'FD'}
            Zii = Xiiminus1_bar;
    end
    U_1Tiiminus1 = U_1Tiiminus1*U{ii};
    Xiiminus1 = reshape(Zii,U_1Tiiminus1,S_iiTd,n);
end
C = reshape(Xiiminus1,U_1Tiiminus1*n,1);
end