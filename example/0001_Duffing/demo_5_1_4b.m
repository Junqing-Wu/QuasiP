% DEMO_5_1_4B  Compute Tor22 using the current three-segment workflow.
% Article: Section 5.1, Fig. 9; P=(2,2), still two base frequencies.
% Segments a and b reuse an integration-based initial guess at Om1=1.94.
% The final bb segment starts at Om1=2.75 using saved qp515_tor22.
% Run demo_5_1_5 before executing this script in full.
% Outputs: qp514_tor22_a, qp514_tor22_b, qp514_tor22_bb.
% The period factors permit doubled tori but do not exclude the Tor11 subspace.
% Next: demo_5_1_4c and figure_5_1_4.
% Setup: execute COCO startup and QuasiP install; work in this example folder.
% COCO outputs normally go to data/<run_name>/. See README.md for the running order.
% Model: qddot + 2*zeta*qdot + q + gamma*q^3
%        = f1*sin(Om1*t) + f2*sin(Om2*t); theta_i=Om_i*t.
% --------------------------------------------------------------------------
rho = [1;1/2^0.5];
Om1 = 1.94;
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
% isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);
isfwithu2 = false;[t2,x2] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[]         ,isfwithu2), tspan, x0, opts);
% isfwithu3 = false; [t3,x3] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing3_F,@Duffing3_L,isfwithu3), tspan, x0, opts);

tt = t2; xx = x2;
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
Xk_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5],[2,2]);
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

%% ====================== tor 22 ======================
prob1 = coco_prob();
% ODE settings: n is the number of model coordinates; d/e count response/excitation frequencies.
% isfwithu indicates whether the nonlinear callback explicitly depends on velocity.
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-4,'h_max',0.1);
% Amplitude settings control selected outputs, velocity reconstruction, and time-domain sampling.
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
% Stability settings: S samples transverse phases, L1 resolves each trajectory,
% N_Ly sets QR iteration count, and Omt selects integration direction(s).
% Near-zero exponents require convergence checks; the stored binary flag alone is insufficient.
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',1000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
% Discretization: MP/mVCF approximates the full torus; MT/SHT integrates trajectory families.
% type selects HB/FD/CO or HT/MF; k sets basis order/stencil; P sets period factors.
% S is the solution sampling grid; S_wave is the denser reconstruction grid.
% OmOd sets MP frequency ordering; Omt and M select SHT integration direction and segment count.
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',    'mVCF'         , ...  
    'isvar'    ,true          , ...
    'type',  { {'HB','HB'} }, ...  
    'k',       { {7,7} }      , ...        
    'P',       { {2,2} }      , ...    
    'S',       { {2^5,2^5} }  , ...
    'S_wave',  { {2^8,2^8} } );
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

% Model callbacks supply nonlinear force, state/parameter derivatives, and the linear part.
% Forced shooting additionally requires derivatives with respect to excitation phases.
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob1, 'qp514_tor22_a', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 1.94]);


prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.1);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',1000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',    'mVCF'         , ...  
    'isvar'    ,true          , ...
    'type',  { {'HB','HB'} }, ...  
    'k',       { {7,7} }      , ...        
    'P',       { {2,2} }      , ...    
    'S',       { {2^5,2^5} }  , ...
    'S_wave',  { {2^8,2^8} } );
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob1, 'qp514_tor22_b', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.94 2.75]);


rho = [1;1/2^0.5];
Om1 = 2.75;%2.71----2.75----2.85
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% Prepare the initial torus guess (integration or saved Fourier data below).
% Direct integration is used only when present; otherwise read a saved torus.
% Read prerequisite saved branch data; the corresponding producer demo must run first.
bd_tor22   = coco_bd_read('qp515_tor22');
MF_f       = coco_bd_col(bd_tor22,'MF_f');
MF_x       = coco_bd_col(bd_tor22,'MF_x');
Xk_MF2_tor22.k_wave = MF_f(:,3:4);
Xk_MF2_tor22.Xc = MF_x;
%----------------------------------------------------------
%% ====================== tor22 COCO ======================
prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2);
prob1 = coco_set(prob1,'cont','PtMX',20);
prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-4,'h_max',0.2);
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {2,2} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^8,2^8} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor22, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob1, 'qp514_tor22_bb', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [2.75 2.8309]);


% Zero constraint preserving the prescribed excitation-frequency ratios.
% Return its analytic Jacobian with respect to the supplied frequency vector.
function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end
