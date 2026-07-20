function Psi = MP_PostP_sta_Psi(data,Om,p,Amp)
n        = data.ode_opts.n;
d        = data.ode_opts.d;
ODE      = data.ODE;
P        = cell2mat(data.PHI.P);
Omi      = data.stable_opts.Omi;
NumOmi   = numel(Omi);
Psi      = cell(2,NumOmi);
if ODE == 2
    N = 2*n;
else
    N = n;
end
isfwithu = data.isfwithu;

for iii = 1:NumOmi
    Om_t     = Omi(iii);
    Om_not   = [Om_t+1:data.ode_opts.d,1:Om_t-1];
    pS       = data.stable_opts.pS{iii};
    if d == 2
        pS = [pS,1];
    end
    pS_prod  = prod(pS);
    pL1      = data.stable_opts.pL1{iii};
    L        = data.stable_opts.L;
    vp       = repmat(p,1,pL1);
    dt1      = 2*pi/(L*Om(Om_t));
    N_Ly     = data.stable_opts.N_Ly;
    %--------------------------------------------------------------------------
    varphi = data.stable_opts.TAU{iii};
    rho    = Om/Om(Om_t);
    if d > 1
        varphi(:,Om_not) = varphi(:,Om_not) + varphi(:,Om_t)*rho(Om_not)';
    end
    %--------------------------------------------------------------------------
    M = cell(pS_prod,1);
    parfor ss = 1:pS_prod
        psii  = (ss-1)*pL1+1:ss*pL1;
        theta = varphi(psii,:)';
        GAM0  = exp(1i*theta'*Amp.MF_f(:,end-d+1:end)');
        X0    = transpose( real(GAM0*Amp.MF_x) );
        if isfwithu
            OMGAM1  = GAM0*1i*diag(Om'*Amp.MF_f(:,end-d+1:end)');
            dX0     = transpose( real(OMGAM1*Amp.MF_x) );
            [DFDX,DFDDX] = data.fx(X0,dX0,vp,theta);
            if isempty(DFDDX)
                DFDDX = zeros(n,n,pL1);
            end
        else
            DFDX  = data.fx(X0,vp,theta);
            DFDDX = zeros(n,n,pL1);
        end
        if ~data.isL
            J = DFDX;
        else
            [~,~,A,B,C] = data.l(p);
            assert(rank(A) == size(A,1), 'Error: Matrix A is singular and cannot be inverted.');
            if ~data.isABC(2)
                B = zeros(n,n);
            end
            if ~data.isABC(3)
                C = zeros(n,n);
            end
            if data.isABC(1)
                J = zeros(2*n,2*n,pL1);
                J(1:n,n+1:2*n,:) = eye(n).*ones(n,n,pL1);
                K = DFDX  - C;
                D = DFDDX - B;
                invM = inv(full(A));
                iMK = pagemtimes(invM,K);
                iMD = pagemtimes(invM,D);
                J(n+1:2*n,1:n,:) = iMK;
                J(n+1:2*n,n+1:2*n,:) = iMD;
            else
                assert(rank(B) == size(B,1), 'Error: Matrix B is singular and cannot be inverted.');
                if ~data.isABC(3)
                    C = zeros(n,n);
                end
                invB = inv(full(B));
                J    = pagemtimes(invB,DFDX-C);
            end
        end
        Mss = J;
        Mss(:,:,1) = eye (N);
        for ii = 1 : pL1-1
            Mss(:,:,ii+1) = expm( (J(:,:,ii)+J(:,:,ii+1))/2*dt1 )*Mss(:,:,ii);
        end
        M{ss} = Mss;
    end
    if d == 1
        Psi = M{1};
    else
        HARM     = data.stable_opts.Harm{iii};
        Mdd = cell2mat(cellfun(@(x) reshape(x(:,:,end),1,[]), M, 'UniformOutput',false));
        B = [];
        M_fun = zeros(pS_prod,N^2);
        for ii = 1:N^2
            Mii = reshape(Mdd(:,ii),pS);
            M_funii = fftn(Mii)/pS_prod;
            M_fun(:,ii) = reshape(M_funii,pS_prod,1);
            max_M_funii = max(abs(M_fun(:,ii)));
            A = find( abs(M_fun(:,ii)) >= max_M_funii*10^(-8) );
            B = union(A, B);
        end
        M_fun  = M_fun(B,:);
        HARM_B = HARM(:,B);

        t = (0:N_Ly-1)'*2*P(Om_t)*pi/Om(Om_t);

        varphi0_ll = zeros(N_Ly,d-1);
        for ii = 1:d-1
            varphi0_ll(:,ii) = mod(t*Om(Om_not(ii)),2*P(Om_not(ii))*pi);
        end
        M_H = exp(1i*varphi0_ll*HARM_B);
        Psi{1,iii} = real(M_H*M_fun);
        Psi{2,iii} = Mdd;        
    end
end
end