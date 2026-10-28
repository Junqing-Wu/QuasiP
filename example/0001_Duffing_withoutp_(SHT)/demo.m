% ddx + 2*zeta*dx + x + gama*x^3 = f1*sin(Om1*t) + f1*sin(Om2*t)
% A(p)*\ddot{x} + D*\dot{x} + K*x = f(x,\dot{x},p,theta) 
%---------------------------------------------------------- 
% duffing1: A = 1; B = 2*zeta; C = 1; 
%           f = -gama*x^3 + f1*sin(theta1) + f1*sin(theta2)
%----------------------------------------------------------
% duffing2: A = []; B = [1 0;0 1]; C = []; 
%           f = [dx
%                -2*zeta*dx - x -gama*x^3 + f1*sin(theta1) + f1*sin(theta2)]
%----------------------------------------------------------
% duffing3: A = []; B = [2*zeta 1;1 0]; C = [1 0;0 -1]; 
%           f = [-gama*x^3 + f1*sin(theta1) + f1*sin(theta2)
%                0]
%----------------------------------------------------------
clearvars;
close all;
clc;
% system parameters
rho = [1;1/2^0.5];
Om1 = 6;
OM0 = rho*Om1;

num_p = 3000;
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
% isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);
isfwithu2 = false;[t2,x2] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@Duffing2_F,[]         ,isfwithu2), tspan, x0, opts);
% isfwithu3 = false; [t3,x3] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@Duffing3_F,@Duffing3_L,isfwithu3), tspan, x0, opts);

tt = t2; xx = x2;
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
ff_final = ff(loc_final,:)/(OM0(1));
X_final  = X(loc_final,:);
X_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5],[2,2]);
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
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',500,'Omi',[1]);
par_opts    = struct('parallel',true,'ncores',12);
prob1 = coco_set(prob1,'cont','PtMX',num_p);
prob1 = coco_set(prob1,'cont','NSV',10);
prob1 = coco_set(prob1,'cont','h_min',0.25,'h_max',0.5);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts,'par_opts',par_opts); 
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @Newmark, ...
    'Omt',        1, ...
    'M_SHT',      1, ...
    'Nsteps',     2^7,...
    'k',          { {5,5} }, ...        
    'P',          { {1,1} }, ...    
    'S',          { {2^5,2^5} }, ...
    'S_wave',     { {2^8,2^8} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing1_F, @Duffing1_DFDX, [], @Duffing1_L, @Duffing1_DFDTH};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd1 = coco(prob1, 'qp1', [], 1, {'k_OM1' 'k_OM2'}, [1.8 6]);
close(11111)

% bd1 = coco_bd_read('qp1');
max_x1    = coco_bd_col(bd1, 'Max_x');
Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1);

figure(1)
plot(Om_mode1(:,1), max_x1,'k--');
hold on
plot(Om_mode1(A11,1), max_x1(A11), 'ko');
xlabel('\omega_1');
ylabel('max(x)');
grid on;

%% ====================== 任务 2 ======================
prob2 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',500,'Omi',[1]);
par_opts    = struct('parallel',true,'ncores',12);
prob2 = coco_set(prob2,'cont','PtMX',num_p);
prob2 = coco_set(prob2,'cont','h_min',0.25,'h_max',0.5);
prob2 = coco_set(prob2,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts,'par_opts',par_opts); 
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @ode45, ...
    'Omt',        1, ...
    'M_SHT',      1, ...
    'Nsteps',     2^8,...
    'k',          { {5,5} }, ...        
    'P',          { {1,1} }, ...    
    'S',          { {2^5,2^5} }, ...
    'S_wave',     { {2^8,2^8} });
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing2_F, @Duffing2_DFDX, [],[], @Duffing2_DFDTH};
prob2 = ode_isol2QuasiP(prob2, '', funcs{:}, X_MF1, OM0);
[data, uidx] = coco_get_func_data(prob2, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob2 = coco_add_func(prob2, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2 = coco(prob2, 'qp2', [], 1, { 'k_OM1' 'k_OM2'}, [1.8 6]);
close(11111)

% bd2 = coco_bd_read('qp2');
max_x2    = coco_bd_col(bd2, 'Max_x');
Om_mode2  = coco_bd_col(bd2, 'Om');
Stab2      = coco_bd_col(bd2, 'Stab');
A2 = find(~Stab2);

figure(1)
plot(Om_mode2(:,1), max_x2(:,1), 'r--');
hold on
plot(Om_mode2(A2,1), max_x2(A2,1), 'r*');
xlabel('\omega_1');
ylabel('max(x)');
grid on;

%% ====================== 任务 3 ======================
prob3 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',500,'Omi',[1]);
par_opts    = struct('parallel',true,'ncores',12);
prob3 = coco_set(prob3,'cont','PtMX',num_p);
prob3 = coco_set(prob3,'cont','h_min',0.25,'h_max',0.5);
prob3 = coco_set(prob3,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts,'par_opts',par_opts); 
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'HT' , ... % HT or MF
    'ODEsolver',  @ode45, ...
    'Omt',        1, ...
    'M_SHT',      1, ...
    'Nsteps',     2^8,...
    'k',          { {5,5} }, ...        
    'P',          { {1,1} }, ...    
    'S',          { {2^5,2^5} }, ...
    'S_wave',     { {2^8,2^8} });
prob3 = coco_set(prob3, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing3_F, @Duffing3_DFDX, [],@Duffing3_L,@Duffing3_DFDTH};
prob3 = ode_isol2QuasiP(prob3, '', funcs{:}, X_MF1, OM0);
[data, uidx] = coco_get_func_data(prob3, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob3 = coco_add_func(prob3, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd3 = coco(prob3, 'qp3', [], 1, { 'k_OM1' 'k_OM2'}, [1.8 6]);
close(11111)

% bd3 = coco_bd_read('qp3');
max_x3    = coco_bd_col(bd3, 'Max_x');
Om_mode3  = coco_bd_col(bd3, 'Om');
Stab3      = coco_bd_col(bd3, 'Stab');
A3 = find(~Stab3);

figure(1)
plot(Om_mode3(:,1), max_x3(:,1), 'b--');
hold on
plot(Om_mode3(A3,1), max_x3(A3,1), 'b^');
xlabel('\omega_1');
ylabel('max(x)');
legend('Addx+Bdx+Cx=f(x,p,t)','unstable','dx=f(x,p,t)','unstable','Bdx+Cx=f(x,p,t)','unstable')
grid on;


function [data, y, J] = fc_2(prob, data, u)
RHO_2 = 1/2^0.5;
Om1 = u(1);
Omi = u(2);
y = RHO_2* Om1 - Omi;
J = [RHO_2 , -1];
end


