% 2ND ODEs: M*\ddot{x} + D*\dot{x} + K*x + n(x,\dot{x}) = e  
% 1ST ODEs: B*\dot{X} = A*X + N(X) + E  ----> X = [] 
% B = [D M; M 0]; A = [-K 0; 0 M]; N = [-n(X); 0]; E = [e; 0];


clearvars;
close all;
clc;
%%
%% ---------------- Parameters (from tables) ----------------
% (Keep both original units and SI where convenient)
Cs10 = 1.08e-10;     % [F]  RF sheath S1 equivalent capacitance constant term
Cs20 = 3.24e-10;     % [F]  Grounded sheath S2 equivalent capacitance constant term

Ii1  = 0.0021;       % [A]  Ion current in sheath S1
Ii2  = 0.0064;       % [A]  Ion current in sheath S2

Ie10 = 0.91;         % [A]  Electron current constant term in sheath S1
Ie20 = 2.74;         % [A]  Electron current constant term in sheath S2

Rpl  = 15.51;        % [Ohm] Plasma bulk equivalent resistance
Lpl  = 2.59e-7;      % [H]   Plasma bulk equivalent inductance

% ---- Geometry / circuit / operating conditions ----
AE_cm2 = 100;        % [cm^2] Powered electrode area
AG_cm2 = 300;        % [cm^2] Ground electrode area
AE     = AE_cm2 * 1e-4;   % [m^2]
AG     = AG_cm2 * 1e-4;   % [m^2]

lbulk_cm = 5.7;      % [cm]   Plasma bulk length
lbulk    = lbulk_cm * 1e-2;  % [m]

Rrf   = 50;          % [Ohm] Equivalent source internal resistance
Rm    = 0.5;         % [Ohm] Equivalent series resistance
Rstray= 0.5;         % [Ohm] Stray parallel branch resistance
Cstray= 200e-12;     % [F]   Stray parallel branch capacitance (200 pF)

Lm2   = 1500e-9;     % [H]   Series matching inductance (1500 nH)
% Cm1   = 1550e-12;    % [F]   Shunt matching capacitance (1550 pF)
% Cm2   = 175e-12;     % [F]   Series matching capacitance (175 pF)

p_mTorr = 5;         % [mTorr] Gas pressure ~= 0.667Pa
p       = 0.667;
Tg      = 300;       % [K]     Gas temperature
Te      = 4.75;      % (as given) Electron temperature
ne      = 1.25e15;   % [m^-3]  Electron density (4.5 W)

Vrf   = 1;         % [V]   RF source voltage amplitude (25 W)
frf   = 13.56e6;     % [Hz]  RF source frequency
wrf   = 2*pi*frf;    % [rad/s]

Cm1   = 1550e-12;    % [F]   Shunt matching capacitance (1550 pF)
Cm2   = 175e-12;     % [F]   Series matching capacitance (175 pF)
%%

% B = CCP_B(p0);
% F0 = CCP_F(tau0,x0,p0,'','');
% 
% X_MF0 = Linear_solution(B,F0.A,F0.E_MF,OM0);


par = [Cm1; Cm2; wrf];x0 = zeros(7,1);
tspan = [0, 1000*2*pi/wrf];

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~,x] = ode15s(@(t,x) ccp_ode_KVL(t,x,par), tspan, x0, opts);

N = 2^10;
dt  = 2*pi/wrf/N;
ccc = 2^8;
[t,x] = ode15s(@(t,x) ccp_ode_KVL(t,x,par), 0:dt:ccc*2*pi/wrf-dt, x(end,:), opts);
bbb = 8;

x = x(end-2^bbb*N+1:end,:);

plot(x(:,1),x(:,4))

X = fft(x)/length(x(:,1));
X(2:end,:)=X(2:end,:)*2;
na = length(x(:,1));
ff = (0:na-1)'/na*2*pi/dt;

figure(4)
plot(ff./wrf , log10(abs(X(:,1))));
xlim([0,15])

XH = ff(1:2^bbb:end,:)./wrf;
XA = X(1:2^bbb:end,:);

%%

OM0 = wrf;
pname = {'cm1' 'cm2'};

p0 = [Cm1;Cm2];
tau0 = zeros(1,1);
X_MF0 = [(0:1:6-1)',XA(1:6,:)];


prob = coco_prob();
ode_opts = struct('d',1, 'e',1 );
prob = coco_set(prob, 'cont','PtMX',1000);
prob = coco_set(prob, 'corr', 'tol', 1e-4);
prob = coco_set(prob, 'QuasiP', 'order', 1, 'DOF', 7, 'ode_opts', ode_opts); 


PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ... 
    'd',       ode_opts.d,...
    'e',       ode_opts.e,...    
    'method',  { {'HB'}   }, ...   
    'isCO1',   0,...                     
    'k',       { {1:1:15} }, ...        
    'U',       { {15*2+1} }, ...    
    'S',       { {2^7}     }, ...
    'S_wave',  { {2^8}     });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 



% PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
%     'type',    'mVCF'           , ... 
%     'd',       ode_opts.d,...
%     'e',       ode_opts.e,...    
%     'method',  { {'HB','FD'}   }, ...   
%     'isCO1',   0,...                     
%     'k',       { {1:1:5,-2:1:1} }, ...        
%     'U',       { {5*2+1,2^5} }, ...    
%     'S',       { {2^5,2^5}     }, ...
%     'S_wave',  { {2^8,2^8}     });
% prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 


% PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
%     'type',    'mVCF'           , ... 
%     'd',       ode_opts.d,...
%     'e',       ode_opts.e,...    
%     'method',  { {'HB',' CO'}   }, ...   
%     'isCO1',   1,...                     
%     'k',       { {1:1:5,[4,1]} }, ...        
%     'U',       { {5*2+1,2^5} }, ...    
%     'S',       { {2^5,2^5}     }, ...
%     'S_wave',  { {2^8,2^8}     });
% prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 


% par_opts.parallel = true;
% par_opts.ncores = 6;
% prob = coco_set(prob, 'QuasiP', 'par_opts', par_opts);

funcs = {@CCP_B, @CCP_B_P, @CCP_F, @CCP_F_X, @CCP_F_P};

prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF0, OM0, pname, p0);

[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');

prob = coco_add_pars(prob, 'OM1', uidx(data.OM_idx(1)), 'Om1');

coco(prob, 'qp', [], 1, { 'Om1' 'cm2' 'cm1'}, [wrf 2*wrf]);







