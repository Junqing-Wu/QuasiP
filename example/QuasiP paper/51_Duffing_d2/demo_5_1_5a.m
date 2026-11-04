clearvars;
close all;
clc;
%% ====================== tor11 ======================
% system parameters
rho = [1;1/2^0.5];
Om1 = 2.85;%2.71----2.75----2.85
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% initial condition
% ---- time intergration method
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
[t2_tor11,x2_tor11] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[] ,isfwithu2), tspan, x0, opts);

tt_tor11 = t2_tor11; xx_tor11 = x2_tor11;
t_tor11 = Num_period/bbb*T1:dt:Num_period*T1-dt;
x_tor11 = zeros(length(t_tor11), size(x0,1));
for j = 1:size(x0,1)
    x_tor11(:,j) = interp1(tt_tor11, xx_tor11(:,j), t_tor11, 'spline');
end
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
X_MF1_tor11 = solve_HOrder(X_final_tor11,ff_final_tor11,rho,[5,5],[1,1]);
X_MF2_tor11 = X_MF1_tor11(:,1:end-1);
%----------------------------------------------------------
% figure(1)
% plot(ff_tor11/(OM0(1)),log(abs(X_tor11(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = X_MF1_tor11(:,1:2)';
% Q   = X_MF1_tor11(:,1+2:end);
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
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2,'Omt',1);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
sta_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',2000,'Omt',[1,2]);
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
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2_tor11, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob1, 'qp515_tor11', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)
%==========================================================================
%% ====================== tor22 ======================
% system parameters
rho = [1;1/2^0.5];
Om1 = 2.75;%2.71----2.75----2.85
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% initial condition
% ---- time intergration method
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
X_MF1_tor22 = solve_HOrder(X_final_tor22,ff_final_tor22,rho,[5,5],[2,2]);
X_MF2_tor22 = X_MF1_tor22(:,1:end-1);
%----------------------------------------------------------
% figure(1)
% plot(ff_tor22/(OM0(1)),log(abs(X_tor22(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = X_MF1_tor22(:,1:2)';
% Q   = X_MF1_tor22(:,1+2:end);
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
sta_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',2000,'Omt',[1,2]);
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
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2_tor22, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob1, 'qp515_tor22', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)
%==========================================================================
%% ====================== tor44 ======================
% system parameters
rho = [1;1/2^0.5];
Om1 = 2.71;%2.71----2.75----2.85
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% initial condition
% ---- time intergration method
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
X_MF1_tor44 = solve_HOrder(X_final_tor44,ff_final_tor44,rho,[5,5],[4,4]);
X_MF2_tor44 = X_MF1_tor44(:,1:end-1);
%----------------------------------------------------------
% figure(1)
% plot(ff_tor44/(OM0(1)),log(abs(X_tor44(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,4])
% k_h = X_MF1_tor44(:,1:2)';
% Q   = X_MF1_tor44(:,1+2:end);
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
sta_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',2000,'Omt',[1,2]);
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
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob1, 'qp515_tor44', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)


function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end
