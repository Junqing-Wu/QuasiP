% DEMO_5_2_P  Compute the periodic branch of the autonomous airfoil/store model.
% Reference: QuasiP article, Section 5.2, Eqs. (40)-(44), Figs. 12-13.
% Coordinates: q1 = relative plunge, q2 = wing pitch, q3 = store pitch.
% Model: Mq*qddot + Cq*qdot + Kq(Q)*q + [0;0;10*q3^3] = 0.
% Integrate at Q=6, estimate the unknown frequency and Fourier coefficients,
% then initialize a d=1 periodic orbit with ode_isol2QuasiP.
% The system is autonomous (e=0); one phase condition fixes the time origin.
% Continue in Q over [5,10] and save qp_d1, including detected TR points.
% The article reports Hopf onset near Q=4.7440, but this script starts
% from integration at Q=6 rather than continuing equilibria from that Hopf point.
% Run this script before the quasi-periodic demos.
% --------------------------------------------------------------------------
clearvars;
close all;
clc;
% system parameters
pname = {'Q'};
p0 = 6;
%%
% initial condition
% ---- time intergration method
x0 = zeros(6,1)+0.1;
Num_period = 2^12;
T1 = 2*pi/1;
S_period = 2^10;
dt = T1/S_period;
bbb = 8;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false; [t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,[],@Airfoil1_F,@Airfoil1_L,isfwithu1), tspan, x0, opts);

tt = t1; xx = x1;
t = Num_period/bbb*T1:dt:Num_period*T1-dt;
x = zeros(length(t), size(x0,1)); 
for j = 1:size(x0,1)
    x(:,j) = interp1(tt, xx(:,j), t, 'spline'); 
end
% Estimate the dominant frequencies and Fourier coefficients from the retained trajectory.
X = fft(x)/length(x(:,1)); X(2:end,:)=X(2:end,:)*2;
na = length(x(:,1)); ff = (0:na-1)'/na*2*pi/dt;
[pks,locs] = findpeaks(log(abs(X(1:end/2,1))));
[pks_sorted, idx_sorted] = sort(pks, 'descend');
keep_idx = pks_sorted >= (pks_sorted(1) - 8);  
pks_final = pks_sorted(keep_idx);
loc_final = locs(idx_sorted(keep_idx));
ff_final = ff(loc_final,:)/1;
X_final  = X(loc_final,:);

% figure(1)
% plot(ff/1,log(abs(X(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,5])

rho = [ff_final(1)];

% Associate spectral peaks with harmonic indices to form the initial Xk_MF.
Xk_MF1 = solve_HOrder(X_final,ff_final,rho,5,1);
Xk_MF2 = Xk_MF1;
Xk_MF2.Xc = Xk_MF1.Xc(:,1:(4-size(Xk_MF1.k_wave,2)));

OM0 = rho;

% % ---- analytical method
% k_h = Xk_MF1.k_wave(:,1)';
% Q   = Xk_MF1.Xc;
% theta = (OM0*t)';
% H_IDF = exp(1i*theta*k_h);
% xx = real(H_IDF*Q);
% 
% figure(2)
% plot(x(:,2),x(:,5))
% hold on
% plot(xx(:,2),xx(:,5),'--')
% xlabel('q2')
% ylabel('dq2')
% legend('TI','TI-FFT')

prob = coco_prob();
% ODE configuration: n=3 physical coordinates; e=0 means autonomous.
% d counts response frequencies, including unknown intrinsic frequencies.
ode_opts = struct('n',3, 'd',1, 'e',0, 'isfwithu',false);
prob = coco_set(prob,'cont','PtMX',1000);
% Reconstruct displacement/velocity output for the selected coordinates.
% Nperiod/Nsampling/Omt determine the requested time-domain sampling.
amp_opts = struct('isamp',true,'isuamp',true,'outdof',[1,2,3],'Nperiod',10,'Nsampling',2^8);
% For d=1 evaluate Floquet stability; for d=2 evaluate Lyapunov stability.
% S/L1/N_Ly control sampling, integration resolution, and iteration length.
sta_opts = struct('issta',true,'L1',2^8);
prob = coco_set(prob,'cont','h_min',1e-3,'h_max',0.05);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts); 
% Choose torus parameterization, numerical method, basis orders, and sampling grids.
% MP/mVCF uses HB/FD/CO; MT/SHT uses HT/MF and the chosen ODE integrator.
% P stores period factors; S_wave controls denser output reconstruction.
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'           , ...     
    'type',    { {'HB'} }, ...   
    'k',       { {9} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^5} }, ...
    'S_wave',  { {2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Airfoil1_F, @Airfoil1_DFDX, [], @Airfoil1_L}; 
% Convert the integration-based Fourier guess to a registered periodic-orbit problem.
prob = ode_isol2QuasiP(prob, '', funcs{:}, Xk_MF2, OM0, pname, p0);
% Continue the one-dimensional branch in Q; inspect actual coverage and stopping reasons.
coco(prob, 'qp_d1', [], 1, { 'Q' 'u_om1'}, [5 10]);