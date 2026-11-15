function [Phi1,Psi,A]  = Newmark(data,ys,phi,phi1s,phi1e,dphi1,Om,p)
%NEWMARK  Integrate one second-order FSE-SHT trajectory segment.
%
% Advance the displacement and velocity along the scaled coordinate phi1
% using the Newmark method. The routine also propagates sensitivities with
% respect to the initial state, base frequencies, and model parameters.
% The terminal state and sensitivity blocks are consumed by MT_SHT_QuasiP.
%
% The physical second-order model is
%   A(p)*ddx + B(p)*dx + C(p)*x = f(x,dx,theta,p).
% With phi1=Om1*t, the scaled equation is
%   Om1^2*A*q'' + Om1*B*q' + C*q = f(q,Om1*q',theta,p).
%
% On input:
%   data   : QuasiP model, PHI settings, and derivative callbacks.
%   ys     : Initial state [q;u] at the beginning of the segment.
%   phi    : Remaining hyper-time phases.
%   phi1s  : Segment start in phi1.
%   phi1e  : Segment end in phi1.
%   dphi1  : Fixed Newmark step in phi1.
%   Om     : Complete base-frequency vector.
%   p      : Model-parameter vector.
%
% On output:
%   Phi1   : Integration grid in phi1.
%   Psi    : State and sensitivity history.
%   A      : Acceleration history.
PHI     = data.PHI;
n       = data.ode_opts.n;
N       = PHI.N;
d       = data.ode_opts.d;
e       = data.ode_opts.e;
a       = data.pdim;
%-------------------------------------------------------------------------------------
Phi1    = (phi1s:dphi1:phi1e)';
n_steps = length(Phi1);
compact = PHI.isvar && isfield(PHI,'endpoint_sensitivities') && PHI.endpoint_sensitivities;
if PHI.isvar && ~compact
    Psi = zeros(n_steps,N+N^2+N*d+N*a);
else
    Psi = zeros(n_steps,N);
end
A = zeros(n_steps,N/2);
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
vs      = us*Om1;
%-------------------------------------------------------------------------------------
Cons    = MT_SHT_get_MDK_Cons(data,p,n,a,Om1);
thetas        = phi;
thetas(Omt)   = phi1s;
thetas(Omnot) = thetas(Omnot) + Rho*phi1s;
if data.isfwithu
    fs  = data.f(qs,vs,p,thetas);
else
    fs  = data.f(qs,p,thetas);
end
as = Cons.Aw\(fs-Cons.Bw*us-Cons.C*qs);
A(1,:) = as';
if PHI.isvar
    Jqs = [eye(n),zeros(n,n),zeros(n,d),zeros(n,a)];
    Jus = [zeros(n,n),eye(n),zeros(n,d),zeros(n,a)];

    [Ms,Ds,Ks,Fs] = MT_SHT_MDKF_var(data,ys,phi1s,phi,Om1,Om2d,Omt,Omnot,Rho,p,Cons,n,d,e,a);

    Jas = Ms\(Fs-Ks*Jqs-Ds*Jus);
end
if PHI.isvar
    Psis     = [qs,Jqs;
                us,Jus];
    if compact
        Psi(1,:) = ys';
    else
        Psi(1,:) = Psis(:);
    end
else
    Psi(1,:) = ys';
end
K_hat = 1/(beta*dphi1^2)*Cons.Aw + gamma/(beta*dphi1)*Cons.Bw +  Cons.C;
for i = 1:n_steps-1

    qe    = qs + us*dphi1 + dphi1^2*(0.5-beta)*as;
    qstar = qs + dphi1*us +(0.5-beta)*dphi1^2*as;
    ustar = us + (1-gamma)*dphi1*as;
    ue    = ustar + gamma/(beta*dphi1)*(qe-qstar);
    ae    = (qe-qstar)/(beta*dphi1^2);
    ve    = Om1*ue;

    phi1  = Phi1(i+1);
    thetae        = phi;
    thetae(Omt)   = phi1;
    thetae(Omnot) = thetae(Omnot) + Rho*phi1;

    ik = 1;
    while 1

        if data.isfwithu
            fe  = data.f(qe,ve,p,thetae);
            [fqe,fve] = data.fx(qe,ve,p,thetae);
            if isempty(fve)
                fve = sparse(n,n);
            end
        else
            fe  = data.f(qe,p,thetae);
            fqe = data.fx(qe,p,thetae);
            fve = sparse(n,n);
        end

        fue = Om1*fve;
        R   = Cons.Aw*ae + Cons.Bw*ue + Cons.C*qe - fe;
        KT  = K_hat - fqe - fue*gamma/(beta*dphi1);
        dqe = -KT\R;

        qe  = qe+dqe;
        ue    = ustar+gamma/(beta*dphi1)*(qe-qstar);
        ae    = (qe-qstar)/(beta*dphi1^2);
        ve    = ue*Om1;
        
        if norm(R)<=RelTol
            break;
        end
        assert(ik<ItMX,'Newton-Raphson does not coverge');
        ik = ik+1;
    end
    if PHI.isvar
       
       [Me,De,Ke,Fe] = MT_SHT_MDKF_var(data,[qe;ue],phi1,phi,Om1,Om2d,Omt,Omnot,Rho,p,Cons,n,d,e,a); 

       Jqstar = Jqs + dphi1*Jus + (0.5-beta)*dphi1^2*Jas;
       Justar = Jus + (1-gamma)*dphi1*Jas;
       S   = Me/(beta*dphi1^2) + gamma*De/(beta*dphi1) + Ke;
       Jqe = S\full(Fe + Me*Jqstar/(beta*dphi1^2) + De*(gamma/(beta*dphi1)*Jqstar-Justar));
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
        if compact
            Psi(i+1,:) = [qs;us]';
        else
            Psi(i+1,:) = Psis(:);
        end
    else
        qs  = qe;
        us  = ue;
        as  = ae;
        ys = [qs;us];
        Psi(i+1,:) = ys';
    end
    A(i+1,:) = as';
end

if compact
    % SHT reads the states at every step and sensitivities at the endpoint.
    terminal = reshape(Psis(:,2:end),1,[]);
    j = 1:numel(terminal);
    Psi = [sparse(Psi),sparse(repmat(n_steps,size(j)),j,terminal,n_steps,numel(j))];
end

% plot(Psi(:,1),Psi(:,2))
% hold on
% plot(Psi(1,1),Psi(1,2),'rx')
% plot(Psi(end,1),Psi(end,2),'bo')
% xlabel('x')

end
