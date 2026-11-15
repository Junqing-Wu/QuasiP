function [data, y, J] = MP_mVCF_QuasiP(prob, data, u)
%MP_MVCF_QUASIP  Evaluate the multi-step variable-coefficient torus residual.
%
% The routine reconstructs physical samples from the coefficient vector Xd,
% applies the discrete hyper-time derivative matrices, evaluates the nonlinear
% force in the sample domain, transforms it back to coefficient space, and
% assembles the residual and Jacobian.
%
% For a second-order model, the linear part represented in the program is
%   R0 = (A*OOmXI2 + B*OmXI1 + C*XI0)*Xd - Theta*Fd.
% For a first-order model, it is
%   R0 = kron(eye(n),OmXI1)*Xd - Fd.
% Here OmXI1=sum_i Om(i)*PHI.XI1{i}, OOmXI2=sum_i,j Om(i)*Om(j)*PHI.XI2{i,j},
% and Fd is the transformed nonlinear force.
%
% On input:
%   prob : COCO problem structure used by callbacks and post-processing.
%   data : Model, PHI discretization, and option data.
%   u    : Unknown vector [Xd(:); Om(:); p(:)].
%
% On output:
%   data : Updated data structure.
%   y    : Discrete m-VCF residual R0.
%   J    : Jacobian with respect to Xd, Om, and p when PHI.isvar is true.
dFddX0 = [];
PHI    = data.PHI;               % PHI basis, transform, and differentiation data
OmOd   = PHI.OmOd;
d      = data.ode_opts.d;        % Number d of response base frequencies
n      = PHI.n;                  % State dimension n
nL     = PHI.nL;                 % Total coefficient count nL=N*L
a      = data.pdim;              % Number a of model parameters p
%----------------------------------------------------------
U       = cell2mat(PHI.U_omod);  % Number U_i of periodic shape functions
U_prod  = prod(U);               % Product U_prod of shape-function counts
pS      = cell2mat(PHI.pS_omod); % Sample counts pS_i=P_i*S_i
pS_prod = prod(pS);              % Product pS_prod of sample counts
%----------------------------------------------------------
useSparse = isfield(PHI,'sparse_jacobian') && PHI.sparse_jacobian && ~data.isfwithu;
R0   = zeros(nL,1);              % NAE residual R0
if PHI.isvar
    if useSparse
        JXd0 = sparse(nL,nL);
    else
        JXd0 = zeros(nL,nL);
    end         % dR0/dC
    Jom0 = zeros(nL,d);          % dR0/dOm
    Jp0  = zeros(nL,a);          % dR0/dp 
else
    JXd0 = [];
    Jom0 = [];
    Jp0  = [];
end
%==========================================================
Xd   = u(1:nL);             % nL*1
Om   = u(nL+1:nL+d);        % d*1
Om   = Om(OmOd);
p    = u(nL+d+1:end);       % a*1
tau  = PHI.tau;             % d*pS_prod
%----------------------------------------------------------
vp   = repmat(u(nL+d+1:end),1,pS_prod);
                            % a*pS_prod
%----------------------------------------------------------
% Reconstruct the torus samples from C=Xd:
%   X0 = Phi(C), and dX0 = dX/dt when velocity dependence is enabled.
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
% Discrete hyper-time derivatives are assembled as
%   D1(Om)=sum_i Om(i)*PHI.XI1{i},
%   D2(Om)=sum_i,j Om(i)*Om(j)*PHI.XI2{i,j}.
XI0 = PHI.XI0;
OmXI1 = Om(1)*PHI.XI1{1};
for ii = 2:d
    OmXI1 = OmXI1 + Om(ii)*PHI.XI1{ii};
end
if PHI.isvar
    OmXI1_om = PHI.XI1;
end
% Continuous model: A(p)*ddx+B(p)*dx+C(p)*x=f(x,dx,theta_1^e,p).
if data.ODE == 2
    OOmXI2    = sparse(U_prod,U_prod);
    for ii = 1:d
        for jj = 1:d
            OOmXI2 = OOmXI2 + Om(ii)*Om(jj)*PHI.XI2{ii,jj};
        end
    end
    if PHI.isvar
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
end
%----------------------------------------------------------
% Transform sampled nonlinear values F0 to coefficient values Fd=S2U(PHI,F0).
% The assembled NAE is R0=linear_operator(Xd,Om,p)-Fd.
if data.isL
    [~,~,A,B,C,Ap,Bp,Cp] = data.l(p);
    if data.isABC(1) == 1
        AOOmXI2 = kron(A,OOmXI2);
        R0     = R0  + AOOmXI2*Xd;
        if PHI.isvar
            JXd0   = JXd0 + AOOmXI2;
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
    end
    if data.isABC(2) == 1
        BOmXI1 = kron(B,OmXI1);
        R0     = R0  + BOmXI1*Xd;
        if PHI.isvar
            JXd0   = JXd0 + BOmXI1;
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
    end
    if data.isABC(3) == 1
        CXI0 = kron(C,XI0);
        R0     = R0  + CXI0*Xd;
        if PHI.isvar
            JXd0   = JXd0 + CXI0;
            if data.isABCwithp(3) == 1
                for aa = 1:a
                    Cpaa = Cp(:,:,aa);
                    if any(reshape(Cpaa,[],1) ~= 0)
                        Jp0(:,aa) = Jp0(:,aa) + kron(Cpaa,XI0)*Xd;
                    end
                end
            end
        end
    end
else
    IOmXI1 = kron(eye(n),OmXI1);
    R0     = R0  + IOmXI1*Xd;
    if PHI.isvar
        JXd0   = JXd0 + IOmXI1;
        for dd = 1:d
            Jom0(:,dd) = Jom0(:,dd) + kron(eye(n),OmXI1_om{dd})*Xd;
        end
    end
end
%----------------------------------------------------------
if data.isfwithu
    F0 = data.f(X0,dX0,vp,tau);
    if PHI.isvar
        [dFdX0,dFddX0] = data.fx(X0,dX0,vp,tau);
    end
    if data.isfwithp
        dFdP0 = data.fp(X0,dX0,vp,tau);
    end
else
    F0 = data.f(X0,vp,tau);
    if PHI.isvar
        dFdX0 = data.fx(X0,vp,tau);
    end
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
    R0 = R0 - PHI.Theta*Fd;
else
    R0 = R0 - Fd;
end
%----------------------------------------------------------
if PHI.isvar
    if ~useSparse
        dFdXd = zeros(nL,nL);
    end
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
    elseif useSparse
        dFdXd = sparse_force_jacobian(PHI,dFdX0,pS_prod,U_prod,nL);
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
        JXd0 = JXd0 - PHI.Theta*dFdXd;
    else
        JXd0 = JXd0 - dFdXd;
    end
end
%----------------------------------------------------------
if data.isfwithp && PHI.isvar
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
y = R0;
if PHI.isvar
    Jom0 = Jom0(:,PHI.pos);
    J = [JXd0, Jom0, Jp0];
else
    J = [];
end
end

function J = sparse_force_jacobian(PHI,df,pS_prod,U_prod,nL)
% Project only coupled DOF pairs, then assemble the global sparse matrix.
[rowDof,colDof] = find(any(df~=0,3));
numBlocks = numel(rowDof);
rows = cell(numBlocks,1);
cols = cell(numBlocks,1);
vals = cell(numBlocks,1);
for iblock = 1:numBlocks
    ii = rowDof(iblock);
    jj = colDof(iblock);
    block = MP_mVCF_S2U_withX(PHI,df(ii,jj,:),pS_prod);
    [i,j,vals{iblock}] = find(block);
    rows{iblock} = i+(ii-1)*U_prod;
    cols{iblock} = j+(jj-1)*U_prod;
end
J = sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),nL,nL);
end

