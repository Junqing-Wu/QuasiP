function Amp = MP_mVCF_PostP_amp(data,Xd,Om,outdof,isuamp)
%MP_MVCF_POSTP_AMP  Reconstruct an m-VCF torus and compute response amplitudes.
%
% This Processor routine converts the coefficient vector Xd to dense
% hyper-time samples, converts those samples to the multi-frequency domain,
% and evaluates the requested amplitude measures over the selected time
% interval. The resulting data are used by PostP_Amp_Sta_bddat.
%
% On input:
%   data    : QuasiP data structure containing PHI and amp_opts.
%   Xd      : m-VCF coefficient vector C for the current solution.
%   Om      : Base-frequency vector omega_1^d.
%   outdof  : Indices of state components whose amplitudes are requested.
%   isuamp  : Flag selecting displacement-only or displacement/velocity output.
%
% On output:
%   Amp     : Amplitude data for the requested state components.

Noutdof  = length(outdof);
PHI      = data.PHI;
n        = PHI.n;
d        = PHI.d;
pwS_omod = cell2mat(PHI.pwS_omod);
pwS_prod = prod(pwS_omod);
U        = cell2mat(PHI.U_omod);
U_prod   = prod(U);
Xd_n     = reshape(Xd,U_prod,n)';
HARM     = transpose(data.amp_opts.Harm);
%-------------------------------------------------------------------------
OmOd     = PHI.OmOd;
pos      = PHI.pos;
HARM     = HARM(:,pos);
%-------------------------------------------------------------------------
Maxht_x  = zeros(1,Noutdof);
maxx  = zeros(1,Noutdof);
minx  = zeros(1,Noutdof);

HT_x     = zeros(pwS_prod,Noutdof);
if isuamp
    HT_u = zeros(pwS_prod,Noutdof);
end
MF_x     = zeros(pwS_prod,Noutdof);
%-------------------------------------------------------------------------
k_address = [];
if d==1
    pwS_omod = [pwS_omod,1];
end
%-------------------------------------------------------------------------
for kk = 1:Noutdof
    ii = outdof(kk);
    if isuamp
        [HT_x(:,kk),HT_u(:,kk)] = MP_mVCF_U2S_xv(Xd_n(ii,:),Om(OmOd),PHI,U_prod);    
    else
        HT_x(:,kk)              = MP_mVCF_U2S_x(Xd_n(ii,:),PHI,U_prod);      
    end
    Maxht_x(kk) = max(abs(HT_x(:,kk)));
    maxx(kk) = max(HT_x(:,kk));
    minx(kk) = min(HT_x(:,kk));
    
    [MF_x(:,kk),A] = HT2MF(HT_x(:,kk),pwS_omod,pwS_prod);
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
%-------------------------------------------------------------------------
end