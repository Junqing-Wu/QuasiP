% DEMO_5_1_6C  Compute and visualize a four-frequency torus.
% Article: Section 5.1, Fig. 11(c); d=e=4 and Om1=2.85.
% Known frequency ratios: [1; 1/sqrt(2); 1/sqrt(5); 1/sqrt(11)].
% No saved input is required. Output: qp516_e4; plotting is included.
% Dense four-dimensional reconstruction can consume substantial memory.
% The current COCO parameter list reports only the first three frequency names.
% Setup: execute COCO startup and QuasiP install; work in this example folder.
% COCO outputs normally go to data/<run_name>/. See README.md for the running order.
% Model: qddot + 0.2*qdot + q + 0.2*q^3 = 5*sum_i sin(Om_i*t).
% --------------------------------------------------------------------------
clearvars;
close all;
clc;
% Define excitation frequencies and model parameters for this experiment.
Om1 = 2.85;
Om2 = Om1/2^0.5;
Om3 = Om1/5^0.5;
Om4 = Om1/11^0.5;

rho = [1;Om2/Om1;Om3/Om1;Om4/Om1];
OM0 = rho*Om1;
%----------------------------------------------------------
% Prepare the initial torus guess (integration or saved Fourier data below).
% Direct integration is used only when present; otherwise read a saved torus.
x0 = [-5.31039662721847;10.7775296427259];
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
% Integrate the physical state to obtain a trajectory for spectral initialization.
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@Duffing_F4,@Duffing_L,isfwithu1), tspan, x0, opts);

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
Xk_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5,5,5],[1,1,1,1]);
Xk_MF2 = Xk_MF1;
Xk_MF2.Xc = Xk_MF1.Xc(:,1:end-1);

%----------------------------------------------------------
% Post-processing: compare reconstructed spectra with integration, then plot torus sections.
% The reshape sizes and section indices below assume the current sampling configuration.
figure(1)
plot(ff/(OM0(1)),log(abs(X(:,1))))
xlabel('normalized f')
ylabel('log(abs(X))')
xlim([0,4])
k_h = Xk_MF1.k_wave(:,1:4)';
Q = Xk_MF1.Xc;
theta = (OM0*t)';
H_IDF = exp(1i*theta*k_h);
xx = real(H_IDF*Q);
figure(2)
plot(x(:,1),x(:,2))
hold on
plot(xx(:,1),xx(:,2),'--')
xlabel('x')
ylabel('dx')
legend('TI','TI-FFT')
%----------------------------------------------------------

%% ====================== d = 4 ======================
prob1 = coco_prob();
% ODE settings: n is the number of model coordinates; d/e count response/excitation frequencies.
% isfwithu indicates whether the nonlinear callback explicitly depends on velocity.
ode_opts    = struct('n',1, 'd',4, 'e',4, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'corr','ItMX',50);
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
    'type',    { {'HB','HB','HB','HB'} }, ...  
    'k',       { {5,3,3,3} }      , ...        
    'P',       { {1,1,1,1} }      , ...    
    'S',       { {2^5,2^5,2^5,2^5} }  , ...
    'S_wave',  { {2^6,2^6,2^6,2^6} }  , ...
    'OmOd' ,   [4,3,2,1]);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

% Model callbacks supply nonlinear force, state/parameter derivatives, and the linear part.
% Forced shooting additionally requires derivatives with respect to excitation phases.
funcs = {@Duffing_F4, @Duffing_DFDX, [], @Duffing_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2, OM0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_4};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Solve the discretized equations at the initial parameter point (COCO dimension 0).
bd4 = coco(prob1, 'qp516_e4', [], 0, {'k_OM1' 'k_OM2' 'k_OM3'}, [1.8 Om1]);


HT_x  = coco_bd_col(bd4, 'HT_x');
HT_u  = coco_bd_col(bd4, 'HT_u');
MF_x  = coco_bd_col(bd4, 'MF_x');
MF_f  = coco_bd_col(bd4, 'MF_f');
T_x  = coco_bd_col(bd4, 'T_x');
T_u  = coco_bd_col(bd4, 'T_u');

figure(11)
plot(ff/(OM0(1)),20*log10(abs(X(:,1))),'k-')
hold on
plot(MF_f(:,1),20*log10(abs(2*MF_x)),'ro','MarkerSize',8)
xlabel('$f$','Interpreter','latex')
ylabel('$\mathrm{dB}$','Interpreter','latex')
xlim([0,4])
ylim([-125,25])
legend('NI','QuasiP','Location','best')
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)

HT_x_4 = reshape(HT_x,2^6,2^6,2^6,2^6);
HT_u_4 = reshape(HT_u,2^6,2^6,2^6,2^6);

figure(22)
hold on
plot(HT_x,HT_u,'.','Color',[0.5 0.5 0.5],'MarkerSize',0.5)
plot(reshape(HT_x_4(1,:,:,:),[],1),reshape(HT_u_4(1,:,:,:),[],1),'.','Color',[0.4 0.4 0.4],'MarkerSize',1)
plot(reshape(HT_x_4(1,1,:,:),[],1),reshape(HT_u_4(1,1,:,:),[],1),'b.','MarkerSize',4)
plot(reshape(HT_x_4(1,1,1,:),[],1),reshape(HT_u_4(1,1,1,:),[],1),'co','MarkerSize',6)
plot(reshape(HT_x_4(1,1,1,1),[],1),reshape(HT_u_4(1,1,1,1),[],1),'g^','MarkerSize',8)
xlabel('$x$','Interpreter','latex','FontSize',30);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',30);
grid on;
legend({'$x(\tau_1,\tau_2,\tau_3,\tau_4)$','$x(\tau_1,\tau_2,\tau_3,0)$','$x(\tau_1,\tau_2,0,0)$','$x(\tau_1,0,0,0)$','$x(0,0,0,0)$'},...
    'Interpreter','latex','Location','best');
legend boxoff;
plot_setting_paper()
xlim([-10,10])
ylim([-28,28])


% Zero constraint preserving the prescribed excitation-frequency ratios.
% Return its analytic Jacobian with respect to the supplied frequency vector.
function [data, y, J] = fc_4(prob, data, u)
rho = [1;1/2^0.5;1/5^0.5;1/11^0.5];
Om1 = u(1);
Om2 = u(2);
Om3 = u(3);
Om4 = u(4);

y = [rho(2)* Om1 - Om2;
     rho(3)* Om1 - Om3;
     rho(4)* Om1 - Om4];
J = [rho(2), -1,  0,  0;
     rho(3),  0, -1,  0;
     rho(4),  0,  0, -1];
end
