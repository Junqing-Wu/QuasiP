%----------------------------------------------------------
clearvars;
close all;
clc;
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
isfwithu2 = false;[t2,x2] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing2_F,[]         ,isfwithu2), tspan, x0, opts);
% isfwithu3 = false; [t3,x3] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing3_F,@Duffing3_L,isfwithu3), tspan, x0, opts);
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
% plot(x(:,1),x(:,2),'c--')
% hold on
% plot(x(1:S_period:end,1),x(1:S_period:end,2),'ro')
% % plot(x(1:2*S_period:end,1),x(1:2*S_period:end,2),'ro')
% % plot(xx(:,1),xx(:,2),'--')
% xlabel('x')
% ylabel('dx')
% % legend('TI','TI-FFT')
%----------------------------------------------------------

%% ====================== 任务 1 ======================
prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2);
par_opts    = struct('parallel',true,'ncores',12);
prob1 = coco_set(prob1,'cont','PtMX',1000);
% prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.5);
stable_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',2000,'Omi',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts,'par_opts',par_opts); 
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @Newmark, ...
    'isvar'    ,  true, ...
    'Omt',        1, ...
    'M_SHT',      8, ...
    'Nsteps',     2^9,...   
    'k',       { {7,7} }, ...        
    'P',       { {2,2} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^7,2^7} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L, @Duffing1_DFDTH};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd1 = coco(prob1, 'qp22', [], 0, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.927 2.83]);
close(11111)

% bd1 = coco_bd_read('qp22');

HT_x = coco_bd_col(bd1,'HT_x');
HT_u = coco_bd_col(bd1,'HT_u');
T_x = coco_bd_col(bd1,'T_x');
T_u = coco_bd_col(bd1,'T_u');

figure(1)
plot(x(:,1),x(:,2),'k-')
hold on
plot(HT_x,HT_u,'r.')
plot(T_x(1:2^2:end),T_u(1:2^2:end),'bo')
legend('TI','HT','kT_1')


MF_x = coco_bd_col(bd1,'MF_x');
MF_f = coco_bd_col(bd1,'MF_f');
tt = coco_bd_col(bd1,'t');



figure(2)
plot(ff/(OM0(1)),log10(abs(X(:,1))),'k-')
hold on
plot(MF_f(:,1),log10(abs(2*MF_x)),'ro')
xlabel('normalized f')
ylabel('log(abs(X))')
xlim([0,6])


function [data, y, J] = fc_2(prob, data, u)
RHO_2 = 1/2^0.5;
Om1 = u(1);
Om2 = u(2);
y = RHO_2* Om1 - Om2;
J = [RHO_2 , -1];
end

