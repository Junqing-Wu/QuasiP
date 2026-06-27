function [data, y, J] = MP_mVCF_QuasiP(prob, data, u)
dFddX0 = [];
PHI    = data.PHI;          % 形函数及其相关矩阵
d      = data.ode_opts.d;   % QP基频的数量
e      = data.ode_opts.e;   % 激励基频的数量
n      = PHI.n;             % QP维度-自由度
nL     = PHI.nL;            % 与QP相关的待定系数的数量
a      = data.pdim;         % continuation paramaters [p] 的数量
%----------------------------------------------------------
U       = cell2mat(PHI.U);  % 各周期形函数的数量
U_prod  = prod(U);          % 周期形函数的数量的乘积
pS      = cell2mat(PHI.pS); % 各基本周期区间内的采样数
pS_prod = prod(pS);         % 基本周期区间内的采样数的乘积
%----------------------------------------------------------
y0   = zeros(nL,1);         % 非线性代数方程组Rd的残差
Jx0  = zeros(nL,nL);        % Rd关于Xd的导数
Jom0 = zeros(nL,d);         % Rd关于Om的导数
Jp0  = zeros(nL,a);         % Rd关于p 的单数 
%==========================================================
Xd   = u(1:nL);             % nL*1
Om   = u(nL+1:nL+d);        % d*1
p    = u(nL+d+1:end);       % a*1
tau  = PHI.tau;             % d*pS_prod
%----------------------------------------------------------
vp   = repmat(u(nL+d+1:end),1,pS_prod);
                            % a*pS_prod
%----------------------------------------------------------
Xd_n = reshape(Xd,U_prod,n)';
if data.isfwithu
    X0  = zeros(n,pS_prod);
    dX0 = zeros(n,pS_prod);
    for nn = 1:n
        [X0(nn,:),dX0(nn,:)] = MP_mVCF_U2S_01(PHI,Xd_n(nn,:),Om,U_prod);
    end
else
    X0  = zeros(n,pS_prod);
    for nn = 1:n
        X0(nn,:) = MP_mVCF_U2S_0(PHI,Xd_n(nn,:),U_prod);
    end
end
%----------------------------------------------------------
XI0 = PHI.XI0;
OmXI1 = Om(1)*PHI.XI1{1};
for ii = 2:d
    OmXI1 = OmXI1 + Om(ii)*PHI.XI1{ii};
end
OmXI1_om = PHI.XI1;
if data.ODE == 2
    OOmXI2    = sparse(U_prod,U_prod);
    for ii = 1:d
        for jj = 1:d
            OOmXI2 = OOmXI2 + Om(ii)*Om(jj)*PHI.XI2{ii,jj};
        end
    end
    OOmXI2_om = cell(1,d);
    for omii = 1:d
        OOmXI2_om{omii} = sparse(U_prod,U_prod);
        for ii = 1:d
            for jj = 1:d
                if omii == ii && omii == jj
                   OOmXI2_om{omii} = OOmXI2_om{omii} + 2*Om(omii)*PHI.XI2{ii,jj};
                elseif omii == ii && omii ~= jj
                   OOmXI2_om{omii} = OOmXI2_om{omii} + Om(jj)*PHI.XI2{ii,jj};
                elseif omii ~= ii && omii == jj
                   OOmXI2_om{omii} = OOmXI2_om{omii} + Om(ii)*PHI.XI2{ii,jj};
                end
            end
        end
    end
end
%----------------------------------------------------------
if data.isL
    [~,~,A,B,C,Ap,Bp,Cp] = data.l(p);
    if data.isABC(1) == 1
        AOOmXI2 = kron(A,OOmXI2);
        y0     = y0  + AOOmXI2*Xd;
        Jx0    = Jx0 + AOOmXI2;
        for dd = 1:d
            Jom0(:,dd) = Jom0(:,dd) + kron(A,OOmXI2_om{dd})*Xd;
        end
        if data.isABCwithp(1) == 1
            for aa = 1:a
                Apaa = Ap(:,:,aa);
                if any(reshape(Apaa,[],1) ~= 0)
                    Jp0(:,aa) = Jp0(:,aa) + kron(Apaa,OOmXI2)*Xd;
                end
            end
        end
    end
    if data.isABC(2) == 1
        BOmXI1 = kron(B,OmXI1);
        y0     = y0  + BOmXI1*Xd;
        Jx0    = Jx0 + BOmXI1;
        for dd = 1:d
            Jom0(:,dd) = Jom0(:,dd) + kron(B,OmXI1_om{dd})*Xd;
        end
        if data.isABCwithp(2) == 1
            for aa = 1:a
                Bpaa = Bp(:,:,aa);
                if any(reshape(Bpaa,[],1) ~= 0)
                    Jp0(:,aa) = Jp0(:,aa) + kron(Bpaa,OmXI1)*Xd;
                end
            end
        end
    end
    if data.isABC(3) == 1
        CXI0 = kron(C,XI0);
        y0     = y0  + CXI0*Xd;
        Jx0    = Jx0 + CXI0;
        if data.isABCwithp(3) == 1
            for aa = 1:a
                Cpaa = Cp(:,:,aa);
                if any(reshape(Cpaa,[],1) ~= 0)
                    Jp0(:,aa) = Jp0(:,aa) + kron(Cpaa,XI0)*Xd;
                end
            end
        end
    end
else
    IOmXI1 = kron(eye(n),OmXI1);
    y0     = y0  + IOmXI1*Xd;
    Jx0    = Jx0 + IOmXI1;
    for dd = 1:d
        Jom0(:,dd) = Jom0(:,dd) + kron(eye(n),OmXI1_om{dd})*Xd;
    end
end
%----------------------------------------------------------
if data.isfwithu
    F0 = data.f(X0,dX0,vp,tau);
    [dFdX0,dFddX0] = data.fx(X0,dX0,vp,tau);
    if data.isfwithp
        dFdP0 = data.fp(X0,dX0,vp,tau);
    end
else
    F0 = data.f(X0,vp,tau);
    dFdX0 = data.fx(X0,vp,tau);
    if data.isfwithp
        dFdP0 = data.fp(X0,vp,tau);
    end
end
%----------------------------------------------------------
Fd = zeros(U_prod,n);
non_zero_nn = find(~all(F0 == 0, 2)); 
for nn = 1:length(non_zero_nn)
    l = non_zero_nn(nn);
    Fd(:,l) = MP_mVCF_S2U(PHI,F0(l,:),pS_prod);
end
Fd = reshape(Fd,[],1);
if PHI.isCO
    y0 = y0 - PHI.Theta*Fd;
else
    y0 = y0 - Fd;
end
%----------------------------------------------------------
dFdXd = zeros(nL,nL);
if data.isfwithu && ~isempty(dFddX0) 
    is_nonzero       = ~all(cat(3, dFdX0,dFddX0) == 0, 3);
    [fx_ll , fx_mm ] = find(is_nonzero);
    for idx = 1:length(fx_ll)
        ll = fx_ll(idx);
        mm = fx_mm(idx);
        if any(dFddX0(ll,mm,:) ~= 0)
            dFdXd_lm = MP_mVCF_S2U_withXDX(PHI, dFdX0(ll,mm,:), ...
                       dFddX0(ll,mm,:), Om, pS_prod);
        else
            dFdXd_lm = MP_mVCF_S2U_withX(PHI, dFdX0(ll,mm,:), pS_prod);
        end
        row_range = (ll-1)*U_prod + 1 : ll*U_prod;
        col_range = (mm-1)*U_prod + 1 : mm*U_prod;
        dFdXd(row_range, col_range) = dFdXd(row_range, col_range) + dFdXd_lm;
    end
%     is_nonzero       = ~all(dFddX0 == 0, 3);
%     dFdOm  = zeros(U_prod,n,d);
%     ddX0dOm = zeros(pS_prod,n,d);
%     [fdx_ll , fdx_mm ] = find(is_nonzero);
%     mm_unique = unique(fdx_mm);
%     for mmu = 1:length(mm_unique)
%         for dd = 1:d
%             ddX0dOm(:,mm_unique(mmu),dd) = MP_mVCF_U2S_ddX0dOm(PHI,Xd_n(mm_unique(mmu),:),U_prod,dd);
%         end
%     end
%     for idx = 1:length(fdx_ll)
%         ll = fdx_ll(idx);
%         mm = fdx_mm(idx);
%         ddX0mmdOm = reshape(ddX0dOm(:,mm,:),[],d);
%         dF0lldOm = diag(reshape(dFddX0(ll,mm,:),[],1))*ddX0mmdOm;
%         for dd = 1:d
%             dFdlldOmdd = MP_mVCF_S2U(PHI,dF0lldOm(:,dd)',pS_prod);
%             dFdOm(:,ll,dd) = dFdOm(:,ll,dd) + dFdlldOmdd;
%         end
%     end
%     dFdOm = reshape(dFdOm,[],d);
%     if PHI.isCO
%         Jom0 = Jom0 - PHI.Theta*dFdOm;
%     else
%         Jom0 = Jom0 - dFdOm;
%     end
else
    is_nonzero     = ~all(dFdX0  == 0, 3);
    [fx_ll, fx_mm] = find(is_nonzero);
    for idx = 1:length(fx_ll)
        ll = fx_ll(idx);
        mm = fx_mm(idx);
        dFdXd_lm = MP_mVCF_S2U_withX(PHI, dFdX0(ll,mm,:), pS_prod);
        row_range = (ll-1)*U_prod + 1 : ll*U_prod;
        col_range = (mm-1)*U_prod + 1 : mm*U_prod;
        dFdXd(row_range, col_range) = dFdXd(row_range, col_range) + dFdXd_lm;
    end
end
if PHI.isCO
    Jx0 = Jx0 - PHI.Theta*dFdXd;
else
    Jx0 = Jx0 - dFdXd;
end
%----------------------------------------------------------
if data.isfwithp
    for aa = 1:a
        Jpaa = zeros(U_prod,n);
        fpaa = dFdP0(:,aa,:);
        fpaa =  reshape(fpaa,n,pS_prod);
        non_zero_nn = find(~all(fpaa == 0, 2));
        for nn = 1:length(non_zero_nn)
            l = non_zero_nn(nn);
            Jpaa(:,l) = MP_mVCF_S2U(PHI,fpaa(l,:),pS_prod);
        end
        Jpaa = reshape(Jpaa,[],1);
        if PHI.isCO
            Jp0(:,aa) = Jp0(:,aa) - PHI.Theta*Jpaa;
        else
            Jp0(:,aa) = Jp0(:,aa) - Jpaa;
        end
    end
end
%==========================================================
y = y0;
J = [Jx0, Jom0, Jp0];
end

