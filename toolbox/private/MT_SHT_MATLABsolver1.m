function [Phi1,Psi,dy]  = MT_SHT_MATLABsolver1(data,ys,phis,phi1s,phi1e,Om,p)
PHI   = data.PHI;
N     = PHI.N;
a     = data.pdim;
d     = data.ode_opts.d;
e     = data.ode_opts.e;
Om1   = Om(PHI.Omt);
Om2d  = Om(PHI.Omnot);
Rho   = Om2d./Om1;

if PHI.isvar
    Psi0  = [ys(:),eye(N),zeros(N,d+a)];
    G          = @(phi1,Psi) MT_SHT_dPSI_function1(data,phi1,Psi,phis,Om1,Om2d,PHI.Omt,PHI.Omnot,Rho,p,N,a,d,e);
else
    Psi0  = ys(:);
    G          = @(phi1,Psi) MT_SHT_dPSI_function2(data,phi1,Psi,phis,Om1,PHI.Omt,PHI.Omnot,Rho,p,N);
end
[Phi1,Psi] = data.PHI.ODEsolver(G, [phi1s phi1e], Psi0(:), data.PHI.ODEsolver_opt);

if data.amp_opts.isamp
    G_fun = @(phi1,Psi) MT_SHT_dPSI_function2(data,phi1,Psi,phis,Om1,PHI.Omt,PHI.Omnot,Rho,p,N);
    nPt   = length(Phi1);
    dy = zeros(nPt, N);
    for k = 1:nPt
        dPsi    = G_fun(Phi1(k), Psi(k,1:N)');
        dy(k,:) = dPsi(:)'*Om1;
    end
else
    dy = [];
end

% plot(Psi(:,1),Psi(:,2))
% hold on
% plot(Psi(1,1),Psi(1,2),'rx')
% plot(Psi(end,1),Psi(end,2),'bo')
% xlabel('x')


end

function dPsi = MT_SHT_dPSI_function1(data,phi1,Psi,phis,Om1,Om2d,Omt,Omnot,Rho,p,N,a,d,e)
%------------------------------------------------------
Psi  = reshape(Psi,N,[]);
y       = Psi(:,1);   
Psi_ys  = Psi(:,2:N+1);
Psi_Om  = Psi(:,N+1+1:N+1+d);
Psi_p   = Psi(:,N+1+d+1:end);
%------------------------------------------------------
phis(Omt)    = phi1;
phis(Omnot)  = phis(Omnot) + phi1*Rho';
theta        = phis(:);
%------------------------------------------------------
if data.ODE == 1
    x  = y;
else
    x  = y(1:N/2);
    dx = y(1+N/2:end);
end
%------------------------------------------------------
if ~data.isL
    g        = data.f(x,p,theta);
    dgdy     = data.fx(x,p,theta);
    if e ~= 0
        dgdth    = data.fth(x,p,theta);
    else
        dgdth    = zeros(N,d);
    end
    if data.isfwithp
       dgdp  = data.fp(x,p,theta);
    end
else
    [~,~,A,B,C,Ap,Bp,Cp] = data.l(p);
    if data.isABC(1) ~= 1
        f        = data.f(x,p,theta);
        dfdx     = data.fx(x,p,theta);
        if e ~= 0
            dfdth    = data.fth(x,p,theta);
        else
            dfdth    = zeros(N,d);
        end
        if data.isABC(3) ~= 1
            C = zeros(N,N);
        end
        g     = B\(f-C*x);
        dgdy  = B\(dfdx-C);
        dgdth = B\dfdth;
        if data.isp
            if data.isABCwithp(2) == 1
                invBp = pagemtimes(pagemtimes(-inv(B),Bp),inv(B));
            else
                invBp = zeros(N,N,a);
            end
            if data.isABCwithp(3) ~= 1
                Cp    = zeros(N,N,a);
            end
            if data.isfwithp
                dfdp  = data.fp(x,p,theta);
            else
                dfdp  = zeros(N,a);
            end
            dgdp  =  reshape(pagemtimes(invBp,f-C*x),[],a) + B\(dfdp - reshape(pagemtimes(Cp,x),[],a) );
        end
    else
        if ~data.isfwithu
            f        = data.f(x,p,theta);
            dfdx     = data.fx(x,p,theta);
            dfddx    = zeros(N/2,N/2);
            if e ~= 0
                dfdth    = data.fth(x,p,theta);
            else
                dfdth    = zeros(N/2,d);
            end
        else
            f        = data.f(x,dx,p,theta);
            [dfdx,dfddx]     = data.fx(x,dx,p,theta);
            if isempty(dfddx)
                dfddx = zeros(N/2,N/2);
            end
            if e ~= 0
                dfdth    = data.fth(x,dx,p,theta);
            else
                dfdth    = zeros(N/2,d);
            end
        end
        if data.isABC(2) ~= 1
            B = zeros(N/2,N/2);
        end
        if data.isABC(3) ~= 1
            C = zeros(N/2,N/2);
        end
        g     = [dx;A\(f-C*x-B*dx)];
        dgdy  = [zeros(N/2,N/2),eye(N/2);
                 A\(dfdx-C)    ,A\(dfddx-B)];
        dgdth = [zeros(N/2,d)  ;A\dfdth];
        if data.isp
            if data.isABCwithp(1) == 1
                invAp = pagemtimes(pagemtimes(-inv(A),Ap),inv(A));
            else
                invAp = zeros(N/2,N/2,a);
            end
            if data.isABCwithp(2) ~= 1
                Bp    = zeros(N/2,N/2,a);
            end
            if data.isABCwithp(3) ~= 1
                Cp    = zeros(N/2,N/2,a);
            end
            if data.isfwithp
                if ~data.isfwithu
                    dfdp  = data.fp(x,p,theta);
                else
                    dfdp  = data.fp(x,dx,p,theta);
                end
            else
                dfdp  = zeros(N/2,a);
            end
            dgdp  =  [zeros(N/2,a);
                     reshape(pagemtimes(invAp,f-C*x-B*dx),[],a) + ...
                     A\(dfdp - reshape(pagemtimes(Cp,x ) ,[],a) + ...
                             - reshape(pagemtimes(Bp,dx) ,[],a) )];
        end
    end
end

dy       = g/Om1;
dPsi_ys  = 1/Om1*dgdy*Psi_ys;
dPsi_Om  = zeros(N,d);
dPsi_Om(:,Omt) = dPsi_Om(:,Omt) + 1/Om1*dgdy*Psi_Om(:,Omt) - 1/Om1^2*g;
for ii = 1:d-1
    dPsi_Om(:,Omt)       = dPsi_Om(:,Omt) - Om2d(ii)/Om1^3*dgdth(:,Omnot(ii))*phi1;
    dPsi_Om(:,Omnot(ii)) = 1/Om1*dgdy*Psi_Om(:,Omnot(ii)) + phi1/Om1^2*dgdth(:,Omnot(ii));
end
if ~data.isp
    dPsi = [dy;dPsi_ys(:);dPsi_Om(:)];
else
    dPsi_p = 1/Om1*dgdy*Psi_p + 1/Om1*dgdp;
    dPsi = [dy;dPsi_ys(:);dPsi_Om(:);dPsi_p(:)];
end

end  



function dPsi = MT_SHT_dPSI_function2(data,phi1,Psi,phis,Om1,Omt,Omnot,Rho,p,N)
%------------------------------------------------------
y       = Psi;
%------------------------------------------------------
phis(Omt)    = phi1;
phis(Omnot)  = phis(Omnot) + phi1*Rho';
theta        = phis(:);
%------------------------------------------------------
if data.ODE == 1
    x  = y;
else
    x  = y(1:N/2);
    dx = y(1+N/2:end);
end
%------------------------------------------------------
if ~data.isL
    g        = data.f(x,p,theta);
else
    [~,~,A,B,C] = data.l(p);
    if data.isABC(1) ~= 1
        f        = data.f(x,p,theta);
        if data.isABC(3) ~= 1
            C = zeros(N,N);
        end
        g     = B\(f-C*x);
    else
        if ~data.isfwithu
            f        = data.f(x,p,theta);
        else
            f        = data.f(x,dx,p,theta);
        end
        if data.isABC(2) ~= 1
            B = zeros(N/2,N/2);
        end
        if data.isABC(3) ~= 1
            C = zeros(N/2,N/2);
        end
        g     = [dx;A\(f-C*x-B*dx)];
    end
end

dPsi  = g/Om1;
end  



