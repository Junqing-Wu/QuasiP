clearvars;
close all;
clc;
% system parameters
rho = [1;1/2^0.5];
Om1 = 5.5;
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% initial condition
x0 = zeros(2,1);
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);

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

%% ====================== 任务 1 ======================
prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',false,'S',2^5,'L',2^8,'N_Ly',500,'Omi',[1,2]);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',   'mVCF'           , ...     
    'type',    { {'HB','HB'} }, ...  
    'k',       { {5,5} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^5,2^5} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd1 = coco(prob1, 'qp1', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [3 5.5]);
close(11111)

% bd1 = coco_bd_read('qp1');
max_x1    = coco_bd_col(bd1, 'Max_x');
Om_mode1  = coco_bd_col(bd1, 'Om');
gamma1    = coco_bd_col(bd1, 'gamma');

figure(1)
plot3(Om_mode1(:,1),gamma1, max_x1,'k--');
ylabel('\gamma_1');
zlabel('\omega_1');
zlabel('max(x)');
grid on;
hold on

%% continuation of QP from previous solution
lab  = coco_bd_labs(bd1, 'EP');
lab  = max(lab);
prob2 = coco_prob();
prob2 = ode_QuasiP2QuasiP(prob2, '', 'qp1', lab);
[data, uidx] = coco_get_func_data(prob2, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob2 = coco_add_func(prob2, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2 = coco(prob2, 'qp2', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)

max_x2    = coco_bd_col(bd2, 'Max_x');
Om_mode2  = coco_bd_col(bd2, 'Om');
gamma2    = coco_bd_col(bd2, 'gamma');

figure(1)
plot3(Om_mode2(:,1),gamma2, max_x2,'k--');
ylabel('\gamma_1');
zlabel('\omega_1');
zlabel('max(x)');
grid on;
hold on

%% continuation of QP from previous solution
lab  = coco_bd_labs(bd2, 'FP');
lab  = max(lab);
prob3 = coco_prob();
prob3 = coco_set(prob3,'cont','PtMX',3000);
prob3 = ode_QuasiP2QuasiP(prob3, '', 'qp2', lab);
[data, uidx] = coco_get_func_data(prob3, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob3 = coco_add_func(prob3, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd3 = coco(prob3, 'qp3', [], 1, {'k_OM1' 'k_OM2' 'gamma' 'zeta' 'f1' 'f2'}, [1.8 5.5]);
close(11111)


max_x3    = coco_bd_col(bd3, 'Max_x');
Om_mode3  = coco_bd_col(bd3, 'Om');
gamma3    = coco_bd_col(bd3, 'gamma');

figure(1)
plot3(Om_mode3(:,1),gamma3, max_x3,'k--');
ylabel('\gamma_1');
zlabel('\omega_1');
zlabel('max(x)');
grid on;
hold on

%% continuation of QP from previous solution
lab  = coco_bd_labs(bd2, 'EP');
lab  = max(lab);
prob4 = coco_prob();
prob4 = coco_set(prob4,'cont','PtMX',3000);
prob4 = ode_QuasiP2QuasiP(prob4, '', 'qp2', lab);
[data, uidx] = coco_get_func_data(prob4, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob4 = coco_add_func(prob4, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd4 = coco(prob4, 'qp4', [], 1, {'k_OM1' 'k_OM2' 'gamma' 'zeta' 'f1' 'f2'}, [1.8 5.5]);
close(11111)

max_x4    = coco_bd_col(bd4, 'Max_x');
Om_mode4  = coco_bd_col(bd4, 'Om');
gamma4    = coco_bd_col(bd4, 'gamma');

figure(1)
plot3(Om_mode4(:,1),gamma4, max_x4,'k--');
ylabel('\gamma_1');
zlabel('\omega_1');
zlabel('max(x)');
grid on;
hold on


function [data, y, J] = fc_2(prob, data, u)
RHO_2 = 1/2^0.5;
Om1 = u(1);
Om2 = u(2);
y = RHO_2* Om1 - Om2;
J = [RHO_2 , -1];
end