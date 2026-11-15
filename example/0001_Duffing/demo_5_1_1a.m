% DEMO_5_1_1A  Compare three equivalent ODE representations.
% Article: Section 5.1, Eqs. (37)-(39), Fig. 6(a).
% Start at Om1=5 and continue in Om1 with Om2=Om1/sqrt(2).
% Build the initial Fourier guess from direct numerical integration.
% No saved run is required. Save Om1=3 events for demo_5_1_1b.
% Outputs: qp511a_sys1, qp511a_sys2, qp511a_sys3.
% Next: demo_5_1_1b, then figure_5_1_1.
% Setup: execute COCO startup and QuasiP install; work in this example folder.
% COCO outputs normally go to data/<run_name>/. See README.md for the running order.
% --------------------------------------------------------------------------
% Two-frequency Duffing oscillator (QuasiP article, Sec. 5.1).
% Physical displacement q(t) satisfies
%   qddot + 2*zeta*qdot + q + gamma*q^3
%       = f1*sin(Om1*t) + f2*sin(Om2*t).
% Here theta1 = Om1*t, theta2 = Om2*t, and p = [zeta; gamma; f1; f2].
% The implemented linear/nonlinear split is
%   A(p)*xddot + B(p)*xdot + C(p)*x = f(x,p,theta).
% -------------------------------------------------------------------------
% Sys. (1), Duffing1: x = q (scalar second-order formulation).
%   A = 1; B = 2*zeta; C = 1;
%   f = -gamma*x^3 + f1*sin(theta1) + f2*sin(theta2).
% -------------------------------------------------------------------------
% Sys. (2), Duffing2: x = [x1; x2] = [q; qdot].
%   A is absent; B = eye(2); C is absent (no L callback is required).
%   xdot = [x2;
%           -2*zeta*x2 - x1 - gamma*x1^3 ...
%           + f1*sin(theta1) + f2*sin(theta2)].
% -------------------------------------------------------------------------
% Sys. (3), Duffing3: x = [x1; x2] = [q; qdot].
%   A is absent; B = [2*zeta, 1; 1, 0]; C = [1, 0; 0, -1];
%   B*xdot + C*x = [-gamma*x1^3 + f1*sin(theta1) + f2*sin(theta2); 0].
%   The second row enforces xdot1 = x2.
% -------------------------------------------------------------------------
clearvars;
close all;
clc;
% Define excitation frequencies and model parameters for this experiment.
Om1 = 5;
Om2 = Om1/2^0.5;
rho = [1;Om2/Om1];
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% Prepare the initial torus guess (integration or saved Fourier data below).
% Direct integration is used only when present; otherwise read a saved torus.
x0 = zeros(2,1);
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
% Integrate the physical state to obtain a trajectory for spectral initialization.
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);
% isfwithu2 = false;[t2,x2] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[]         ,isfwithu2), tspan, x0, opts);
% isfwithu3 = false; [t3,x3] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing3_F,@Duffing3_L,isfwithu3), tspan, x0, opts);

tt = t1; xx = x1;
t = Num_period/bbb*T1:dt:Num_period*T1-dt;
x = zeros(length(t), size(x0,1)); 
for j = 1:size(x0,1)
    x(:,j) = interp1(tt, xx(:,j), t, 'spline'); 
end

% Estimate Fourier coefficients from the retained integration window.
% ff is an angular-frequency grid; this is an initialization estimate, not a residual check.
X = fft(x)/length(x(:,1)); X(2:end,:)=X(2:end,:)*2;
na = length(x(:,1)); ff = (0:na-1)'/na*2*pi/dt;
[pks,locs] = findpeaks(log(abs(X(1:end/2,1))));
[pks_sorted, idx_sorted] = sort(pks, 'descend');
keep_idx = pks_sorted >= (pks_sorted(1) - 8);  
pks_final = pks_sorted(keep_idx);
loc_final = locs(idx_sorted(keep_idx));
loc_final = [1;loc_final];

ff_final = ff(loc_final,:)/(OM0(1));
X_final  = X(loc_final,:);
% Match spectral peaks to integer/fractional combinations of the prescribed base frequencies.
% Xk_MF stores harmonic indices k_wave and their complex amplitudes Xc.
Xk_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5],[1,1]);
Xk_MF2 = Xk_MF1;
Xk_MF2.Xc = Xk_MF1.Xc(:,1:end-1);

%----------------------------------------------------------
% figure(1)
% plot(ff/(OM0(1)),log(abs(X(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = Xk_MF1.k_wave(:,1:2)';
% Q   = Xk_MF1.Xc;
% theta = (OM0*t)';
% H_IDF = exp(1i*theta*k_h);
% xx = real(H_IDF*Q);
% figure(2)
% plot(x(:,1),x(:,2))
% hold on
% plot(xx(:,1),xx(:,2),'--')
% xlabel('x')
% ylabel('dx')
% legend('TI','TI-FFT')
%----------------------------------------------------------

%% ====================== sys 1 ======================
prob1 = coco_prob();
% ODE settings: n is the number of model coordinates; d/e count response/excitation frequencies.
% isfwithu indicates whether the nonlinear callback explicitly depends on velocity.
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);

% Amplitude settings control selected outputs, velocity reconstruction, and time-domain sampling.
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
% Stability settings: S samples transverse phases, L1 resolves each trajectory,
% N_Ly sets QR iteration count, and Omt selects integration direction(s).
% Near-zero exponents require convergence checks; the stored binary flag alone is insufficient.
sta_opts = struct('issta',false,'S',2^5,'L1',2^8,'N_Ly',500,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
% Discretization: MP/mVCF approximates the full torus; MT/SHT integrates trajectory families.
% type selects HB/FD/CO or HT/MF; k sets basis order/stencil; P sets period factors.
% S is the solution sampling grid; S_wave is the denser reconstruction grid.
% OmOd sets MP frequency ordering; Omt and M select SHT integration direction and segment count.
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',    'mVCF'         , ...  
    'type',    { {'HB','HB'} }, ...  
    'k',       { {5,5} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^5,2^5} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

% Model callbacks supply nonlinear force, state/parameter derivatives, and the linear part.
% Forced shooting additionally requires derivatives with respect to excitation phases.
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
prob1 = coco_add_event(prob1, 'Om1=3', 'k_OM1', 3);


% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob1, 'qp511a_sys1', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.8 Om1]);
close(11111)
%% ====================== sys 2 ======================
prob2 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
sta_opts = struct('issta',false,'S',2^6,'L1',2^8,'N_Ly',500,'Omt',[1,2]);
prob2 = coco_set(prob2,'cont','PtMX',3000);
prob2 = coco_set(prob2,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'           , ...     
    'type',    { {'HB','HB'} }, ...   
    'k',       { {5,5} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} }, ...
    'OmOd' ,   [1,2]);
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing2_F, @Duffing2_DFDX, @Duffing2_DFDP};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob2 = ode_isol2QuasiP(prob2, '', funcs{:}, Xk_MF1, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob2, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob2 = coco_add_func(prob2, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
prob2 = coco_add_event(prob2, 'Om1=3', 'k_OM1', 3);
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob2, 'qp511a_sys2', [], 1, { 'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.8 Om1]);
close(11111)
%% ====================== sys 3 ======================
prob3 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
sta_opts = struct('issta',false,'S',2^6,'L1',2^8,'N_Ly',500,'Omt',[1,2]);
prob3 = coco_set(prob3,'cont','PtMX',3000);
prob3 = coco_set(prob3,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',   'mVCF'           , ...     
    'type',    { {'HB','HB'} }, ...   
    'k',       { {5,5} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} }, ...
    'OmOd' ,   [1,2]);
prob3 = coco_set(prob3, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing3_F, @Duffing3_DFDX, @Duffing3_DFDP,@Duffing3_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob3 = ode_isol2QuasiP(prob3, '', funcs{:}, Xk_MF1, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob3, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob3 = coco_add_func(prob3, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
prob3 = coco_add_event(prob3, 'Om1=3', 'k_OM1', 3);
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob3, 'qp511a_sys3', [], 1, { 'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.8 Om1]);
close(11111)

% Zero constraint preserving the prescribed excitation-frequency ratios.
% Return its analytic Jacobian with respect to the supplied frequency vector.
function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end






