clearvars;
close all;
clc;
% system parameters
Om1 = 2.85;
Om2 = Om1/2^0.5;
rho = [1;Om2/Om1];
OM0 = rho*Om1;
%----------------------------------------------------------
% initial condition
% ---- time intergration method
x0 = [-5.31039662721847;10.7775296427259];
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@Duffing_F2,@Duffing_L,isfwithu1), tspan, x0, opts);

tt = t1; xx = x1;
t = Num_period/bbb*T1:dt:Num_period*T1-dt;
x = zeros(length(t), size(x0,1)); 
for j = 1:size(x0,1)
    x(:,j) = interp1(tt, xx(:,j), t, 'spline'); 
end

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
X_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5],[1,1]);
X_MF2 = X_MF1(:,1:end-1);

%----------------------------------------------------------
% figure(1)
% plot(ff/(OM0(1)),log(abs(X(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = X_MF1(:,1:2)';
% Q   = X_MF1(:,1+2:end);
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

%% ====================== 任务 1 ======================
prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'corr','ItMX',50);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',false,'S',2^5,'L',2^8,'N_Ly',500,'Omi',[1,2]);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',    'mVCF'         , ...  
    'type',    { {'HB','HB'} }, ...  
    'k',       { {7,5} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^5,2^5} }  , ...
    'S_wave',  { {2^6,2^6} }  , ...
    'OmOd' ,   [2,1]);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@Duffing_F2, @Duffing_DFDX, [], @Duffing_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2 = coco(prob1, 'qp516_e3', [], 0, {'k_OM1' 'k_OM2' 'k_OM3'}, [1.8 Om1]);

HT_x  = coco_bd_col(bd2, 'HT_x');
HT_u  = coco_bd_col(bd2, 'HT_u');
MF_x  = coco_bd_col(bd2, 'MF_x');
MF_f  = coco_bd_col(bd2, 'MF_f');
T_x  = coco_bd_col(bd2, 'T_x');
T_u  = coco_bd_col(bd2, 'T_u');

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

HT_x_3 = reshape(HT_x,2^6,2^6);
HT_u_3 = reshape(HT_u,2^6,2^6);

figure(22)
hold on
plot(reshape(HT_x_3(:,:),[],1),reshape(HT_u_3(:,:),[],1),'b.','MarkerSize',4)
plot(reshape(HT_x_3(1,:),[],1),reshape(HT_u_3(1,:),[],1),'co','MarkerSize',6)
plot(reshape(HT_x_3(1,1),[],1),reshape(HT_u_3(1,1),[],1),'g^','MarkerSize',8)
xlabel('$x$','Interpreter','latex','FontSize',30);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',30);
grid on;
legend({'$x(\tau_1,\tau_2)$','$x(\tau_1,0)$','$x(0,0)$'},...
    'Interpreter','latex','Location','best');
legend boxoff;
plot_setting_paper()
xlim([-10,10])
ylim([-25,25])



function [data, y, J] = fc_2(prob, data, u)
rho = [1;1/2^0.5];
Om1 = u(1);
Om2 = u(2);

y = rho(2)* Om1 - Om2;
J = [rho(2), -1];
end
