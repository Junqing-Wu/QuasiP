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
Om1 = 5;
Om2 = Om1/2^0.5;
rho = [1;Om2/Om1];
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
isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);
% isfwithu2 = false;[t2,x2] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[]         ,isfwithu2), tspan, x0, opts);
% isfwithu3 = false; [t3,x3] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing3_F,@Duffing3_L,isfwithu3), tspan, x0, opts);

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

%% ====================== sys 1 ======================
prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);

amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',false,'S',2^5,'L',2^8,'N_Ly',500,'Omi',[1,2]);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',    'mVCF'         , ...  
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
prob1 = coco_add_event(prob1, 'Om1=3', 'k_OM1', 3);


coco(prob1, 'qp511a_sys1', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.8 Om1]);
close(11111)
%% ====================== sys 2 ======================
prob2 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',false,'S',2^6,'L',2^8,'N_Ly',500,'Omi',[1,2]);
prob2 = coco_set(prob2,'cont','PtMX',3000);
prob2 = coco_set(prob2,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts); 
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
prob2 = ode_isol2QuasiP(prob2, '', funcs{:}, X_MF1, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob2, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob2 = coco_add_func(prob2, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
prob2 = coco_add_event(prob2, 'Om1=3', 'k_OM1', 3);
coco(prob2, 'qp511a_sys2', [], 1, { 'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.8 Om1]);
close(11111)
%% ====================== sys 3 ======================
prob3 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',false,'S',2^6,'L',2^8,'N_Ly',500,'Omi',[1,2]);
prob3 = coco_set(prob3,'cont','PtMX',3000);
prob3 = coco_set(prob3,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts); 
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
prob3 = ode_isol2QuasiP(prob3, '', funcs{:}, X_MF1, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob3, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob3 = coco_add_func(prob3, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
prob3 = coco_add_event(prob3, 'Om1=3', 'k_OM1', 3);
coco(prob3, 'qp511a_sys3', [], 1, { 'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.8 Om1]);
close(11111)

function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end






