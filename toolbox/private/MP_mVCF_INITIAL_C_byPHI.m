function C = MP_mVCF_INITIAL_C_byPHI(data,Xk_MF)
%MP_MVCF_INITIAL_C_BYPHI  Convert multi-frequency initial data to m-VCF coefficients.
%
% Xk_MF stores harmonic indices k_wave and state coefficients Xc. The routine
% first evaluates this representation on the PHI hyper-time samples, then
% successively projects each hyper-time direction into the coefficient ordering
% used by the m-VCF tensor-product basis.
%
% On input:
%   data  : QuasiP data containing PHI and the number of state variables.
%   Xk_MF : Initial multi-frequency representation with fields k_wave and Xc.
%
% On output:
%   C     : Column vector in the coefficient ordering expected by QuasiP_init_sol.
PHI = data.PHI;
n   = data.ode_opts.n;

pS   = PHI.pS_omod;
U    = PHI.U_omod;
d    = PHI.d;
tau  = PHI.tau;

k_h = Xk_MF.k_wave';
Q   = Xk_MF.Xc;

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