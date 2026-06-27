function [Phi1,Psi]  = Newmark(data,ys,phi,phi1s,phi1e,dphi1,Om,p)
% 0: M*a + D*u + K*q = f(q,Om1*u,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
%    M = Om1^2*A; D = Om1*B 
PHI     = data.PHI;
n       = data.ode_opts.n;
N       = PHI.N;
d       = data.ode_opts.d;
e       = data.ode_opts.e;
a       = data.pdim;
%-------------------------------------------------------------------------------------
Phi1    = (phi1s:dphi1:phi1e)';
n_steps = length(Phi1);
if PHI.isvar
    Psi = zeros(n_steps,N+N^2+N*d+N*a);
else
    Psi = zeros(n_steps,N);
end
%-------------------------------------------------------------------------------------
ItMX    = PHI.ODE2slover_opt.ItMX;
RelTol  = PHI.ODE2slover_opt.RelTol;
alpha   = PHI.ODE2slover_opt.alpha;
beta    = (1+alpha)^2/4;
gamma   = 1/2+alpha;
%-------------------------------------------------------------------------------------
Omt     = PHI.Omt;
Omnot   = PHI.Omnot;
Om1     = Om(Omt);
Om2d    = Om(Omnot);
Rho     = Om2d/Om1;
%-------------------------------------------------------------------------------------
qs      = ys(1:n);
us      = ys(n+1:2*n);
vs      = Om1*us;
%-------------------------------------------------------------------------------------
Cons    = MT_SHT_get_MDK_Cons(data,p,Om1,n,a);
thetas        = phi;
thetas(Omt)   = phi1s;
thetas(Omnot) = thetas(Omnot) + phi(Omnot)+Rho*phi1s;
if data.isfwithu
    fs  = data.f(qs,vs,p,thetas);
else
    fs  = data.f(qs,p,thetas);
end
as = Cons.Aw\(fs-Cons.C*qs-Cons.Bw*us);
if PHI.isvar
    Jqs = [eye(n),zeros(n,n),zeros(n,d),zeros(n,a)];
    Jus = [zeros(n,n),eye(n),zeros(n,d),zeros(n,a)];
    [Ms,Ds,Ks,Fs] = MT_SHT_MDKF_var(data,ys,phi1s,phi,Om1,Om2d,Omt,Omnot,Rho,p,Cons,n,d,e,a);
    Jas = Ms\(Fs-Ks*Jqs-Ds*Jus);
end
if PHI.isvar
    Psis     = [qs,Jqs;
                us,Jus];
    Psi(1,:) = Psis(:);
else
    Psi(1,:) = ys';
end

for i = 1:n_steps-1
    phi1  = Phi1(i+1);
    qe    = qs;
    qstar = qs+dphi1*us+(0.5-beta)*dphi1^2*as;
    ustar = us+(1-gamma)*dphi1*as;
    ue    = ustar+gamma/(beta*dphi1)*(qe-qstar);
    ae    = (qe-qstar)/(beta*dphi1^2);
    thetae        = phi;
    thetae(Omt)   = phi1;
    thetae(Omnot) = thetae(Omnot) + phi(Omnot)+Rho*phi1s;
    ik = 1;
    while 1
        if data.isfwithu
            fe  = data.f(qe,ve,p,thetae);
            [fqe,fve] = data.fx(qe,ve,p,thetae);
            if isempty(fve)
                fve = zeros(n,n);
            end
        else
            fe  = data.f(qe,p,thetae);
            fqe = data.fx(qe,p,thetae);
            fve = zeros(n,n);
        end
        fue = Om1*fve;
        r   = Cons.Aw*ae + Cons.Bw*ue + Cons.C*qe - fe;
        S   = Cons.Aw/(beta*dphi1^2) + gamma*(Cons.Bw - fue)/(beta*dphi1) + (Cons.C-fqe);
        dqe = -S\r;
        qe  = qe+dqe;
        ue  = ustar+gamma/(beta*dphi1)*(qe-qstar);
        ae  = (qe-qstar)/(beta*dphi1^2);
        if norm(dqe)<=RelTol*norm(qe)
            break;
        end
        assert(ik<ItMX,'Newton-Raphson does not coverge');
        ik = ik+1;
    end

    if PHI.isvar
       [Me,De,Ke,Fe] = MT_SHT_MDKF_var(data,[qe;ue],phi1,phi,Om1,Om2d,Omt,Omnot,Rho,p,Cons,n,d,e,a); 
       Jqstar = Jqs+dphi1*Jus+(0.5-beta)*dphi1^2*Jas;
       Justar = Jus+(1-gamma)*dphi1*Jas;
       S   = Me/(beta*dphi1^2) + gamma*De/(beta*dphi1) + Ke;
       Jqe = S\(Fe + Me*Jqstar/(beta*dphi1^2) + De*(gamma/(beta*dphi1)*Jqstar-Justar));
       Jue = Justar+gamma/(beta*dphi1)*(Jqe-Jqstar);
       Jae = (Jqe-Jqstar)/(beta*dphi1^2);
    end
    if PHI.isvar
        qs  = qe;
        us  = ue;
        as  = ae;
        Jqs = Jqe;
        Jus = Jue;
        Jas = Jae;
        Psis  = [qs,Jqs;
                 us,Jus];
        Psi(i+1,:) = Psis(:);
    else
        qs  = qe;
        us  = ue;
        as  = ae;
        ys = [qs;us];
        Psi(i+1,:) = ys';
    end
end

% plot(Psi(:,1),Psi(:,4))
% hold on
% plot(Psi(1,1),Psi(1,4),'rx')
% plot(Psi(end,1),Psi(end,4),'bo')

end


