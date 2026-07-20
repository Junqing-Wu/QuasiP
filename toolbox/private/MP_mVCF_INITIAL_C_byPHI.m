function C = MP_mVCF_INITIAL_C_byPHI(data,X_MF)

PHI = data.PHI;
n   = data.ode_opts.n;

pS   = PHI.pS_omod;
U    = PHI.U_omod;
d    = PHI.d;
tau  = PHI.tau;

k_h = X_MF(:,1:d)';
Q   = X_MF(:,1+d:end);

H_IDF = exp(1i*tau'*k_h);
Xiiminus1 = real(H_IDF*Q);

% compute discrete Zii
pS_iiTd = prod(cell2mat(pS));
U_1Tiiminus1 = 1;
Xiiminus1 = reshape(Xiiminus1,U_1Tiiminus1,pS_iiTd,n);

for ii = 1:d
    pS_iiTd = pS_iiTd/pS{ii};
    Xiiminus1_bar = pagetranspose(reshape(Xiiminus1,U_1Tiiminus1,pS{ii},pS_iiTd,n) );
    iGAM = PHI.iGam0{ii};
    switch PHI.me_omod{ii}
        case{'HB'}
            Zii = pagemtimes(iGAM,Xiiminus1_bar);
        case{'CO'}
            Zii = Xiiminus1_bar;           
        case{'FD'}
            Zii = Xiiminus1_bar;
    end
    U_1Tiiminus1 = U_1Tiiminus1*U{ii};
    Xiiminus1 = reshape(Zii,U_1Tiiminus1,pS_iiTd,n);
end
C = reshape(Xiiminus1,U_1Tiiminus1*n,1);
end