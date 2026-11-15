% DEMO_5_1_5  Prepare representative tori, sections, and spectra.
% Article: Section 5.1, Fig. 10. No saved input is required.
% Tor11: Om1=2.85, P=(1,1); Tor22: Om1=2.75, P=(2,2);
% Tor44: Om1=2.71, P=(4,4). Each COCO call is a single-point solve.
% Outputs: qp515_tor11, qp515_tor22, qp515_tor44 and their _b variants.
% The paired runs use different Omt choices for time-domain reconstruction.
% Keep the direct-integration arrays in the workspace and run figure_5_1_5 immediately.
% The base Tor22/Tor44 runs also provide initial guesses for demos 4b/4c.
% Setup: execute COCO startup and QuasiP install; work in this example folder.
% COCO outputs normally go to data/<run_name>/. See README.md for the running order.
% Model: qddot + 2*zeta*qdot + q + gamma*q^3
%        = f1*sin(Om1*t) + f2*sin(Om2*t); theta_i=Om_i*t.
% --------------------------------------------------------------------------
clearvars;
close all;
clc;
%% ====================== tor11 ======================
% Define excitation frequencies and model parameters for this experiment.
rho = [1;1/2^0.5];
Om1 = 2.85;%2.71----2.75----2.85
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
isfwithu2 = false;
[t2_tor11,x2_tor11] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[] ,isfwithu2), tspan, x0, opts);

tt_tor11 = t2_tor11; xx_tor11 = x2_tor11;
t_tor11 = Num_period/bbb*T1:dt:Num_period*T1-dt;
x_tor11 = zeros(length(t_tor11), size(x0,1));
for j = 1:size(x0,1)
    x_tor11(:,j) = interp1(tt_tor11, xx_tor11(:,j), t_tor11, 'spline');
end
% Estimate Fourier coefficients from the retained integration window.
% ff is an angular-frequency grid; this is an initialization estimate, not a residual check.
X_tor11 = fft(x_tor11)/length(x_tor11(:,1)); X_tor11(2:end,:)=X_tor11(2:end,:)*2;
na = length(x_tor11(:,1)); ff_tor11 = (0:na-1)'/na*2*pi/dt;
[pks_tor11,locs_tor11] = findpeaks(log(abs(X_tor11(1:end/2,1))));
[pks_sorted_tor11, idx_sorted_tor11] = sort(pks_tor11, 'descend');
keep_idx_tor11 = pks_sorted_tor11 >= (pks_sorted_tor11(1) - 8);
pks_final_tor11 = pks_sorted_tor11(keep_idx_tor11);
loc_final_tor11 = locs_tor11(idx_sorted_tor11(keep_idx_tor11));
loc_final_tor11 = [1;loc_final_tor11];
ff_final_tor11 = ff_tor11(loc_final_tor11,:)/(OM0(1));
X_final_tor11  = X_tor11(loc_final_tor11,:);
% Match spectral peaks to integer/fractional combinations of the prescribed base frequencies.
% Xk_MF stores harmonic indices k_wave and their complex amplitudes Xc.
Xk_MF1_tor11 = solve_HOrder(X_final_tor11,ff_final_tor11,rho,[5,5],[1,1]);
Xk_MF2_tor11 = Xk_MF1_tor11;
Xk_MF2_tor11.Xc = Xk_MF1_tor11.Xc(:,1:end-1);
%----------------------------------------------------------
% figure(1)
% plot(ff_tor11/(OM0(1)),log(abs(X_tor11(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = Xk_MF1_tor11.k_wave(:,1:2)';
% Q   = Xk_MF1_tor11.Xc;
% theta = (OM0*t_tor11)';
% H_IDF = exp(1i*theta*k_h);
% xx = real(H_IDF*Q);
% figure(2)
% plot(x_tor11(:,1),x_tor11(:,2),'c--')
% hold on
% plot(x_tor11(1:S_period:end,1),x_tor11(1:S_period:end,2),'ro')
% % plot(x(1:2*S_period:end,1),x(1:2*S_period:end,2),'ro')
% % plot(xx(:,1),xx(:,2),'--')
% xlabel('x')
% ylabel('dx')
% % legend('TI','TI-FFT')
%----------------------------------------------------------
%% ====================== tor11 COCO ======================
prob1 = coco_prob();
% ODE settings: n is the number of model coordinates; d/e count response/excitation frequencies.
% isfwithu indicates whether the nonlinear callback explicitly depends on velocity.
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
% Amplitude settings control selected outputs, velocity reconstruction, and time-domain sampling.
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',1);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
% Stability settings: S samples transverse phases, L1 resolves each trajectory,
% N_Ly sets QR iteration count, and Omt selects integration direction(s).
% Near-zero exponents require convergence checks; the stored binary flag alone is insufficient.
sta_opts = struct('issta',true,'S',2^6,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
% Discretization: MP/mVCF approximates the full torus; MT/SHT integrates trajectory families.
% type selects HB/FD/CO or HT/MF; k sets basis order/stencil; P sets period factors.
% S is the solution sampling grid; S_wave is the denser reconstruction grid.
% OmOd sets MP frequency ordering; Omt and M select SHT integration direction and segment count.
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {1,1} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^5,2^5} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
% Model callbacks supply nonlinear force, state/parameter derivatives, and the linear part.
% Forced shooting additionally requires derivatives with respect to excitation phases.
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor11, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
coco(prob1, 'qp515_tor11', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)

prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
sta_opts = struct('issta',true,'S',2^6,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {1,1} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^5,2^5} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor11, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
coco(prob1, 'qp515_tor11_b', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)
%==========================================================================
%% ====================== tor22 ======================
% Define excitation frequencies and model parameters for this experiment.
rho = [1;1/2^0.5];
Om1 = 2.75;%2.71----2.75----2.85
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
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
% isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);
isfwithu2 = false;
[t2_tor22,x2_tor22] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[] ,isfwithu2), tspan, x0, opts);

tt_tor22 = t2_tor22; xx_tor22 = x2_tor22;
t_tor22 = Num_period/bbb*T1:dt:Num_period*T1-dt;
x_tor22 = zeros(length(t_tor22), size(x0,1));
for j = 1:size(x0,1)
    x_tor22(:,j) = interp1(tt_tor22, xx_tor22(:,j), t_tor22, 'spline');
end
X_tor22 = fft(x_tor22)/length(x_tor22(:,1)); X_tor22(2:end,:)=X_tor22(2:end,:)*2;
na = length(x_tor22(:,1)); ff_tor22 = (0:na-1)'/na*2*pi/dt;
[pks_tor22,locs_tor22] = findpeaks(log(abs(X_tor22(1:end/2,1))));
[pks_sorted_tor22, idx_sorted_tor22] = sort(pks_tor22, 'descend');
keep_idx_tor22 = pks_sorted_tor22 >= (pks_sorted_tor22(1) - 8);
pks_final_tor22 = pks_sorted_tor22(keep_idx_tor22);
loc_final_tor22 = locs_tor22(idx_sorted_tor22(keep_idx_tor22));
loc_final_tor22 = [1;loc_final_tor22];
ff_final_tor22 = ff_tor22(loc_final_tor22,:)/(OM0(1));
X_final_tor22  = X_tor22(loc_final_tor22,:);
Xk_MF1_tor22 = solve_HOrder(X_final_tor22,ff_final_tor22,rho,[5,5],[2,2]);
Xk_MF2_tor22 = Xk_MF1_tor22;
Xk_MF2_tor22.Xc = Xk_MF1_tor22.Xc(:,1:end-1);
%----------------------------------------------------------
% figure(1)
% plot(ff_tor22/(OM0(1)),log(abs(X_tor22(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = Xk_MF1_tor22.k_wave(:,1:2)';
% Q   = Xk_MF1_tor22.Xc;
% theta = (OM0*t_tor22)';
% H_IDF = exp(1i*theta*k_h);
% xx = real(H_IDF*Q);
% figure(2)
% plot(x_tor22(:,1),x_tor22(:,2),'c--')
% hold on
% plot(x_tor22(1:S_period:end,1),x_tor22(1:S_period:end,2),'ro')
% % plot(x(1:2*S_period:end,1),x(1:2*S_period:end,2),'ro')
% % plot(xx(:,1),xx(:,2),'--')
% xlabel('x')
% ylabel('dx')
% % legend('TI','TI-FFT')
%----------------------------------------------------------
%% ====================== tor22 COCO ======================
prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',1);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
sta_opts = struct('issta',true,'S',2^6,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {2,2} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^5,2^5} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor22, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
coco(prob1, 'qp515_tor22', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)

prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
sta_opts = struct('issta',true,'S',2^6,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {2,2} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^5,2^5} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor22, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
coco(prob1, 'qp515_tor22_b', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)
%==========================================================================
%% ====================== tor44 ======================
% Define excitation frequencies and model parameters for this experiment.
rho = [1;1/2^0.5];
Om1 = 2.71;%2.71----2.75----2.85
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
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
% isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);
isfwithu2 = false;
[t2_tor44,x2_tor44] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[] ,isfwithu2), tspan, x0, opts);

tt_tor44 = t2_tor44; xx_tor44 = x2_tor44;
t_tor44 = Num_period/bbb*T1:dt:Num_period*T1-dt;
x_tor44 = zeros(length(t_tor44), size(x0,1));
for j = 1:size(x0,1)
    x_tor44(:,j) = interp1(tt_tor44, xx_tor44(:,j), t_tor44, 'spline');
end
X_tor44 = fft(x_tor44)/length(x_tor44(:,1)); X_tor44(2:end,:)=X_tor44(2:end,:)*2;
na = length(x_tor44(:,1)); ff_tor44 = (0:na-1)'/na*2*pi/dt;
[pks_tor44,locs_tor44] = findpeaks(log(abs(X_tor44(1:end/2,1))));
[pks_sorted_tor44, idx_sorted_tor44] = sort(pks_tor44, 'descend');
keep_idx_tor44 = pks_sorted_tor44 >= (pks_sorted_tor44(1) - 8);
pks_final_tor44 = pks_sorted_tor44(keep_idx_tor44);
loc_final_tor44 = locs_tor44(idx_sorted_tor44(keep_idx_tor44));
loc_final_tor44 = [1;loc_final_tor44];
ff_final_tor44 = ff_tor44(loc_final_tor44,:)/(OM0(1));
X_final_tor44  = X_tor44(loc_final_tor44,:);
Xk_MF1_tor44 = solve_HOrder(X_final_tor44,ff_final_tor44,rho,[5,5],[4,4]);
Xk_MF2_tor44 = Xk_MF1_tor44;
Xk_MF2_tor44.Xc = Xk_MF1_tor44.Xc(:,1:end-1);
%----------------------------------------------------------
% figure(1)
% plot(ff_tor44/(OM0(1)),log(abs(X_tor44(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = Xk_MF1_tor44.k_wave(:,1:2)';
% Q   = Xk_MF1_tor44.Xc;
% theta = (OM0*t_tor44)';
% H_IDF = exp(1i*theta*k_h);
% xx = real(H_IDF*Q);
% figure(2)
% plot(x_tor44(:,1),x_tor44(:,2),'c--')
% hold on
% plot(x_tor44(1:S_period:end,1),x_tor44(1:S_period:end,2),'ro')
% % plot(x(1:2*S_period:end,1),x(1:2*S_period:end,2),'ro')
% % plot(xx(:,1),xx(:,2),'--')
% xlabel('x')
% ylabel('dx')
% % legend('TI','TI-FFT')
%----------------------------------------------------------
%% ====================== tor44 COCO ======================
prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',1);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
sta_opts = struct('issta',true,'S',2^6,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {4,4} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^5,2^5} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
coco(prob1, 'qp515_tor44', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)

prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
sta_opts = struct('issta',true,'S',2^6,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {4,4} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^5,2^5} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
coco(prob1, 'qp515_tor44_b', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)

% Zero constraint preserving the prescribed excitation-frequency ratios.
% Return its analytic Jacobian with respect to the supplied frequency vector.
function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end
