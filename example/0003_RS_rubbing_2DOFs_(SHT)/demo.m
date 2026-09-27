% A(p)*\ddot{x} + D*\dot{x} + K*x = f(x,\dot{x},p,theta) 
%---------------------------------------------------------- 
% RS1: A = diag([m,m]); B = diag([cs,cs]); C = diag([ks,ks]); 
%           f1 = - kc*Theta*(1-r0/r)*(X-mu*Nu*Y) + m*e*OM^2*cos(theta1);
%           f2 = - kc*Theta*(1-r0/r)*(mu*Nu*X+Y) + m*e*OM^2*sin(theta1);
%----------------------------------------------------------
% RS2: A = []; B = eye(4); C = []; 
%           f1 = U;
%           f2 = V;
%           f3 = ( -cs*U - ks*X - kc*Theta*(1-r0/r)*(X-mu*Nu*Y) + m*e*OM^2*cos(theta1) )/m;
%           f4 = ( -cs*V - ks*Y - kc*Theta*(1-r0/r)*(mu*Nu*X+Y) + m*e*OM^2*sin(theta1) )/m;
%%
clearvars;
close all;
clc;
% system parameters
rho = 1;
Om1 = 0.3;
OM0 = rho*Om1;
pname = {'mu','OM'};
p0 = [0.2;OM0(1)];
%%
% initial condition
% ---- time intergration method
x0 = zeros(4,1)+0.0001;
Num_period = 2^9;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
% isfwithu1 = true; [t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@RS1_F,@RS1_L,isfwithu1), tspan, x0, opts);
isfwithu2 = false;[t2,x2] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@RS2_F,[]    ,isfwithu2), tspan, x0, opts);

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
 
rho = ff_final(1:2);
X_MF = solve_HOrder(X_final,ff_final,rho,[5,5],[1,1]);
X_MF2 = X_MF(:,1:4);
% ---- analytical method
OM0 = rho*Om1;

% figure(1)
% plot(ff/(OM0(1)),log(abs(X(:,1))))
% xlabel('normalized f')
% ylabel('log(abs(X))')
% xlim([0,15])
% k_h = X_MF(:,1:2)';
% Q   = X_MF(:,1+2:end);
% theta = (OM0*t)';
% H_IDF = exp(1i*theta*k_h);
% xx = real(H_IDF*Q);
% 
% figure(2)
% plot(x(:,1),x(:,2))
% hold on
% plot(xx(:,1),xx(:,2),'--')
% xlabel('x')
% ylabel('dx')
% legend('TI','TI-FFT')


prob1 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',1, 'isfwithu',true);
amp_opts = struct('isamp',true,'outdof',[1 2],'Nperiod',1,'Nsampling',2^8);
par_opts = struct('parallel',true,'ncores',12);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'cont','h_min',.1,'h_max',1);

prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'par_opts',par_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MT', ...
    'type',       'SHT'           , ...     
    'method',     'HT' , ... % HT or MF
    'ODEsolver',  @ode45, ...
    'isvar'    ,  true, ...
    'Omt',        1, ...
    'M_SHT',      4, ...
    'Nsteps',     2^10,...   
    'k',       { {7,7} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^6,2^6} }, ...
    'S_wave',  { {2^8,2^8} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@RS1_F, @RS1_DFDX, @RS1_DFDP, @RS1_L, @RS1_DFDTH};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname , p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
prob1 = coco_add_glue(prob1,'OMp2',uidx(data.kOM_idx),uidx(data.p_idx(2)));
bd1 = coco(prob1, 'qp1', [], 1, { 'k_OM1' 'u_om2' 'OM' 'mu'}, [0.3 1]);
close(11111)

figure(1)
OM1 = coco_bd_col(bd1,'k_OM1');
om2 = coco_bd_col(bd1,'u_om2');
plot(OM1,om2);
xlabel('k_{OM1}');
ylabel('\omega_2');
grid on;

max_x = coco_bd_col(bd1, 'Max_x');
figure(2)
plot(OM1,max_x(:,1));
xlabel('k_{OM1}');
ylabel('max(x)');
grid on;




prob2 = coco_prob();
ode_opts = struct('n',4, 'd',2, 'e',1, 'isfwithu',false);
prob2 = coco_set(prob2,'cont','PtMX',3000);
prob2 = coco_set(prob2,'QuasiP','ode_opts',ode_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB'} }, ...   
    'k',       { {7,7} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^8,2^8} }, ...
    'S_wave',  { {2^8,2^8} });
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@RS2_F, @RS2_DFDX, @RS2_DFDP};
prob2 = ode_isol2QuasiP(prob2, '', funcs{:}, X_MF, OM0, pname , p0);
[data, uidx] = coco_get_func_data(prob2, 'QuasiP', 'data', 'uidx');
prob2 = coco_add_glue(prob2,'OMp2',data.kOM_idx,data.p_idx(2));

bd2 = coco(prob2, 'qp11', [], 1, { 'k_OM1' 'u_om2' 'OM' 'mu'}, [0.3 1]);
close(11111)

figure(1) 
hold on
OM1 = coco_bd_col(bd2,'k_OM1');
om2 = coco_bd_col(bd2,'u_om2');
plot(OM1,om2,'--')      

