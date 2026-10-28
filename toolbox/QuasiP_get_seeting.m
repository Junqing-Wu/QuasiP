function data = QuasiP_get_seeting(prob,tbid,data)

% QUASIP_GET_SETTINGS   Read 'QuasiP' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% DATA = QUASIP_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data strcture.
%--------------------------------------------------------------------------------------------------------------------
ode_opts.n           = [];            % Number of DOFs (must be provided by user)
ode_opts.d           = [];            % Number of base frequencies of QP (must be provided by user)
ode_opts.e           = [];            % Number of base frequencies of Excitation (must be provided by user)
ode_opts.isfwithu    = false;         % Default: f(x,p,theta1,...,thetae); Conversely: f(x,dx,p,theta1,...,thetae)
defaults.ode_opts    = ode_opts;      % Store ODE options in defaults
%--------------------------------------------------------------------------------------------------------------------
par_opts.parallel    = false;         % Enable or disable parallel computing
cl = parcluster('local');
par_opts.ncores      = cl.NumWorkers; % Number of CPU cores for parallel computation
defaults.par_opts    = par_opts;      % Store parallel settings in default configuration
%--------------------------------------------------------------------------------------------------------------------
amp_opts.isamp       = false;         % Enable or disable amplitude calculation
amp_opts.isuamp      = false;       
amp_opts.isht        = false;
amp_opts.outdof      = 1;             % DOFs for output (all or specific indices)
amp_opts.Nperiod     = 1;             % Number of periods for time domain: [0,Nperiod*Ti]
amp_opts.Nsampling   = 2^8;           % Number of sampling points within [0,Ti)
amp_opts.Omt         = 1;             % The i of period Ti
amp_opts.t0_per      = 0;             % t0 = t0_per*Ti
defaults.amp_opts    = amp_opts;      % Store amplitude options in default configuration
%--------------------------------------------------------------------------------------------------------------------
stable_opts.issta    = false;         % Enable or disable stability analysis
stable_opts.Omi      = 1;
stable_opts.S        = 2^5;           % Number of trajectories along a single hyper-time direction
stable_opts.L        = 2^8;           % Number of sampling points for integration
stable_opts.N_Ly     = 1000;          % Iterations for Lyapunov exponents (when d > 1)
defaults.stable_opts = stable_opts;   % Store stability analysis options in default configuration
%--------------------------------------------------------------------------------------------------------------------
copts = coco_get(prob, tbid);
copts = coco_merge(defaults, copts);
%--------------------------------------------------------------------------------------------------------------------
assert(~isempty(copts.ode_opts.n),'%s: The ode_opts is not defined');
data.ode_opts = copts.ode_opts;
%--------------------------------------------------------------------------------------------------------------------
if ~isfield(copts, 'PHI_INITIAL') % defaults.method: VCHB
    types = repmat({'HB'}, 1, data.ode_opts.d);
    k = repmat({5}    , 1, data.ode_opts.d);
    P = repmat({1}    , 1, data.ode_opts.d);
    S = repmat({2^5}  , 1, data.ode_opts.d);
    S_w = repmat({2^8}, 1, data.ode_opts.d);
    OmOd = 1:data.ode_opts.d;
    copts.PHI_INITIAL = struct( 'param', 'MP', ...
        'method',  'mVCF'           , ...
        'type',    { types   }, ...
        'k',       { k    }, ...
        'P',       { P    }, ...
        'S',       { S    }, ...
        'S_wave',  { S_w  },... 
        'OmOd',    OmOd, ...
        'd',       data.ode_opts.d, ...
        'n',       data.ode_opts.n);
else
    copts.PHI_INITIAL.d = data.ode_opts.d;
    copts.PHI_INITIAL.n = data.ode_opts.n;
end
PHI = COM_PHI(copts.PHI_INITIAL);
data.PHI = PHI;
%--------------------------------------------------------------------------------------------------------------------
if copts.amp_opts.isamp
    copts.amp_opts.Harm =  fun_k(cell2mat(PHI.P_omod),cell2mat(PHI.wS_omod));
    if data.ode_opts.d <= 2
        copts.amp_opts.isht = true;
    end
else
    copts.stable_opts.issta = false;
end
data.amp_opts   =  copts.amp_opts;
%--------------------------------------------------------------------------------------------------------------------
if copts.stable_opts.issta
    copts.stable_opts = fill_stable_TAU_data(copts.stable_opts, data.PHI, data.ode_opts);
    if prod(copts.stable_opts.pS{1}) >1 && ~copts.par_opts.parallel && strcmp(PHI.param, 'MP')
       copts.par_opts.parallel = true;
    end
end
data.stable_opts  = copts.stable_opts;
%--------------------------------------------------------------------------------------------------------------------
if copts.par_opts.parallel
    % activate parallel computing
    activate_parallel(copts.par_opts.ncores)
end
data.par_opts   = copts.par_opts;
end
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
function stable_opts = fill_stable_TAU_data(stable_opts, PHI, ode_opts)
S        = stable_opts.S;
L        = stable_opts.L;
period   = cell2mat(PHI.P);
Omi      = stable_opts.Omi;
d        = ode_opts.d;
assert(length(unique(Omi)) == length(Omi) && all(Omi <= d), ...
   '%s: The OmiforStab is not defined');
Num_Omi = numel(Omi);
TAU      = cell(1,Num_Omi);
pL1      = cell(1,Num_Omi);
Harm     = cell(1,Num_Omi);
pS       = cell(1,Num_Omi);
for ii = 1:Num_Omi
    Om_t    = Omi(ii);
    Om_not  = [Om_t+1:d,1:Om_t-1];
    per_t   = period(Om_t);
    per_not = period(Om_not);
    tau_1im1 = 0;
    pS_1im1  = 1;
    for jj = 1:d-1
        tau_i = (0:2*pi/S:per_not(jj)*2*pi-2*pi/S)';
        tau_1im1 = [ kron( ones(per_not(jj)*S,1),tau_1im1 ) , kron( tau_i,ones(pS_1im1,1) ) ];
        pS_1im1 = pS_1im1*per_not(jj)*S;
    end
    tau_t = (0:2*pi/L:per_t*2*pi)';
    TAUii = zeros(pS_1im1*(per_t*L+1),d);
    TAUii(:,Om_t) = kron(ones(pS_1im1,1),tau_t);
    TAUii(:,Om_not) = kron(tau_1im1(:,2:end),ones(per_t*L+1,1));
    TAU{ii} = TAUii;
    pL1{ii} = per_t*L+1;
    if d > 1
        if ~isfield(PHI,'Omt')
            Harm{ii} = fun_k(per_not,repmat(S,1,d-1) );
        else
            if Om_t == PHI.Omt
                my_S = cell2mat(PHI.S);
                Harm{ii} = fun_k(per_not, my_S(Om_not));
            else
                Harm{ii} = fun_k(per_not,repmat(S,1,d-1) );
            end
        end
        pS{ii}   = per_not.*repmat(S,1,d-1);
    else
        Harm{ii} = [];
        pS{ii}   = [];
    end
end
stable_opts.TAU  = TAU;
stable_opts.pL1  = pL1;
stable_opts.Harm = Harm;
stable_opts.pS   = pS;
end
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
function k = fun_k(p,wS)
pwS = p.*wS;
d = length(p);
k = zeros(d,prod(pwS));
for jj = 1:d

    k1 = 0:1/p(jj):wS(jj)-1/p(jj);
    k2 = (0:1/p(jj):wS(jj)-1/p(jj))-wS(jj);
    kk = zeros(1,pwS(jj));
    for ii = 1:pwS(jj)
        if abs(k1(ii))<abs(k2(ii))
            kk(ii) = k1(ii);
        elseif abs(k1(ii))>abs(k2(ii))
            kk(ii) = k2(ii);
        else
            kk(ii) = -abs(k1(ii));
        end
    end
    N1 = prod(pwS(1:jj-1));
    N2 = prod(pwS(jj+1:end));
    k(jj,:) = kron(ones(1,N2),kron(kk,ones(1,N1)));
end
end
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
function activate_parallel(varargin)
%ACTIVATE_PARALLEL This function starts parallel pool for the first time
%and detects the number of cores avaliable
try
    pp1 = gcp('nocreate');
    if isempty(pp1)
        h = helpdlg('Starting parallel pool for the first time and detecting number of available cores.', 'Info');
        disp('----------------------------------------------------------------')
        disp('Starting parallel pool for the first time and detecting number')
        disp('of available cores.')
        disp('----------------------------------------------------------------')
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        if numel(varargin)>=1 && isnumeric(varargin{1})
            parpool(myCluster,varargin{1});
        else
            parpool(myCluster);
        end
        pp2 = gcp('nocreate');
        cpuNum =pp2.NumWorkers;
        save('cluster_info.mat','cpuNum')

        if isvalid(h)
            close(h)
        end
    end
catch ME
    rethrow(ME)
end
end
