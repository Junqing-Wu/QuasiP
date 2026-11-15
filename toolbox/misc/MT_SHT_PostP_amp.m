function Amp = MT_SHT_PostP_amp(data,Om,outdof,isuamp)
%MT_SHT_POSTP_AMP  Compute amplitudes from stored FSE-SHT trajectories.
%
% MT_SHT_QuasiP stores trajectory states and optional velocities when amplitude
% processing is enabled. This routine interpolates the trajectory samples,
% transforms the remaining hyper-time directions, and evaluates amplitudes over
% the time interval selected by amp_opts.
%
% On input:
%   data    : QuasiP data containing PHI, stored trajectory samples, and amp_opts.
%   Om      : Base-frequency vector.
%   outdof  : State-component indices to report.
%   isuamp  : Flag selecting displacement-only or displacement/velocity output.
%
% On output:
%   Amp     : Amplitude data for the requested state components.

Noutdof   = length(outdof);
PHI       = data.PHI;
N         = PHI.N;
n         = PHI.n;
d         = PHI.d;
Omt       = PHI.Omt;
Omnot     = PHI.Omnot;
Om1       = Om(Omt);
Om2d      = Om(Omnot);
Rho       = Om2d./Om1;
pydy_temp = data.phi_y_dy;
pS2d      = PHI.pS2d;
Nsteps    = PHI.L1;
wS1       = PHI.S_wave{Omt}*PHI.P{Omt};
if d>1
    pwS2d     = prod(cell2mat(PHI.pwSnot));
else
    pwS2d     = 1;
end
pwS_omod  = cell2mat(PHI.pwS_omod);
pwS_prod  = prod(pwS_omod);
HARM      = transpose(data.amp_opts.Harm);
wVarphi1  = PHI.wVarphi1;
%-------------------------------------------------------------------------
pos       = PHI.pos;
HARM      = HARM(:,pos);
%-------------------------------------------------------------------------
Maxht_x   = zeros(1,Noutdof);
maxx      = zeros(1,Noutdof);
minx      = zeros(1,Noutdof);
HT_x      = zeros((wS1+1)*pwS2d,Noutdof);
if isuamp
    HT_u  = zeros((wS1+1)*pwS2d,Noutdof);
end
MF_x      = zeros( wS1*pwS2d   ,Noutdof);
%-------------------------------------------------------------------------
if PHI.IsFixedStep
    ydy = pydy_temp(1,2:end,1:pS2d);
    for jj = 1:PHI.M
        st = (jj-1)*pS2d + 1;
        ed = jj*pS2d;
        ydy = [ydy;pydy_temp(2:end,2:end,st:ed)];
    end
    ydy  = PHI.SPLine*reshape(ydy,Nsteps*PHI.P{Omt}+1,[]);
    ydy  = reshape(ydy,wS1+1,N,pS2d);
    ydy(:,n+1:2*n,:) = ydy(:,n+1:2*n,:)*Om1;
else
    ydy = cellfun(@(x) x(1,:), pydy_temp(1:pS2d), 'UniformOutput',false);
    for jj = 1:PHI.M
        st = (jj-1)*pS2d + 1;
        ed = jj*pS2d;
        ydy = cellfun(@(a,b)[a;b(2:end,:)], ydy, pydy_temp(st:ed), 'UniformOutput',false);
    end
    ydy = cellfun(@(dat) interp1(dat(:,1), dat(:,2:end), wVarphi1, 'pchip', 'extrap'),...
        ydy, 'UniformOutput',false);
    ydy = cat(3, ydy{:});
end
%-------------------------------------------------------------------------
k_address = [];
if d==1
    pwS_omod = [pwS_omod,1];
end
%-------------------------------------------------------------------------
if d>1
    wGOiwG = cell(wS1,1);
    kQ     = PHI.kQ;
    wGam   = PHI.wGam;
    iwGam  = PHI.iwGam;
    parfor ss = 2:wS1+1
        lamda = wVarphi1(ss)*Rho'*kQ;
        v0     = [1,kron(cos(lamda),[1,1])];
        v1     = kron(sin(lamda),[0,1]);
        O_phi1 = sparse(diag(v0)+diag(-v1,1)+diag(v1,-1));
        wGOiwG{ss,1} = wGam*O_phi1*iwGam;
    end
end
%-------------------------------------------------------------------------
for kk = 1:Noutdof
    ii = outdof(kk);
    x_ii_temp  = reshape(ydy(:,ii,:),wS1+1,pS2d);
    if d>1
        x_ii       = PHI.wGamiGam*x_ii_temp';
        HT_x(:,kk) = reshape(x_ii',[],1);
        if isuamp
            if data.ODE == 2
                u_ii_temp = reshape(ydy(:,n+ii,:),wS1+1,pS2d);
            else
                u_ii_temp = reshape(ydy(:,N+ii,:),wS1+1,pS2d);
            end
            u_ii      = transpose(PHI.wGamiGam*u_ii_temp');
            HT_u(:,kk) = reshape(u_ii,[],1);
        end
        x_tau = x_ii;
        parfor ss = 2:wS1+1
            x_tau(:,ss) = wGOiwG{ss,1}*x_tau(:,ss);
        end
        x_tau = x_tau';
        x_tau = reshape(x_tau(1:end-1,:),[],1);
    else
        HT_x(:,kk) = x_ii_temp;
        x_tau      = x_ii_temp(1:end-1);
        if isuamp
            if data.ODE == 2
                u_ii_temp = reshape(ydy(:,n+ii,:),wS1+1,pS2d);
            else
                u_ii_temp = reshape(ydy(:,N+ii,:),wS1+1,pS2d);
            end
            HT_u(:,kk) = u_ii_temp;
        end
    end
    Maxht_x(kk) = max(abs(HT_x(:,kk)));
    maxx(kk) = max(HT_x(:,kk));
    minx(kk) = min(HT_x(:,kk));
    [MF_x(:,kk),A] = HT2MF(x_tau,pwS_omod,pwS_prod);
    k_address = union(A, k_address);
end
%-------------------------------------------------------------------------
MF_x = MF_x(k_address,:);
k    = HARM(k_address,:);
kom  = k*Om;
MF_f = [kom/Om(1),kom,k];
%-------------------------------------------------------------------------
[T_x, T_u, tspan] = deal([]);
%-------------------------------------------------------------------------
if d~=1 && data.amp_opts.Nperiod~=0
    T1 = 2*pi/Om(data.amp_opts.Omt);
    dt = T1/data.amp_opts.Nsampling;
    t  = data.amp_opts.t0_per*T1:dt:(data.amp_opts.t0_per+data.amp_opts.Nperiod)*T1;
    theta = (Om*t)';
    qt_GAM = exp(1i*theta*k');
    T_x = real(qt_GAM*MF_x);
    if isuamp
        qt_GAMdot = qt_GAM*1i*diag(kom);
        T_u = real(qt_GAMdot*MF_x);
    end
    tspan = t;
end
if d == 1
    T1 = 2*pi/Om(1);
    Sampling = cell2mat(PHI.S_wave);
    PPP      = cell2mat(PHI.P);
    dt       = T1/Sampling;
    tspan = 0:dt:PPP*T1-dt;
end
%-------------------------------------------------------------------------
Amp.Maxht_x = Maxht_x;
Amp.maxx    = maxx;
Amp.minx    = minx;

Amp.HT_x    = HT_x;
Amp.T_x     = T_x;
if isuamp
    Amp.HT_u = HT_u;
    Amp.T_u  = T_u;
end
Amp.MF_x    = MF_x;
Amp.MF_f    = MF_f;
Amp.t       = tspan;
end