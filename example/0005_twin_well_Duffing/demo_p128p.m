% ddx + 2*zeta*dx + x + gama*x^3 = f1*sin(Om1*t) + f2*sin(Om2*t)
% A(p)*\ddot{x} + D*\dot{x} + K*x = f(x,\dot{x},p,theta) 
%---------------------------------------------------------- 
% duffing1: A = 1; B = 2*zeta; C = 1; 
%           f = -gama*x^3 + f1*sin(theta1) + f2*sin(theta2) + f3*sin(theta3)
%----------------------------------------------------------
% clearvars;
% close all;
% clc;
% system parameters
rho = 1;
Om1 = 2.15467;
OM0 = rho*Om1;
pname = {'r' 'b' 'f'};
p0 = [0.5;1;1];
%----------------------------------------------------------
% initial condition
% ---- time intergration method
x0 = [1;0];
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^9;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false;[t1,x1] = ode45(@(t,x) my_dxdt(t,x,p0,OM0,@Duffing1_F,@Duffing1_L,isfwithu1), tspan, x0, opts);

tt = t1; xx = x1;
t = Num_period/bbb*T1:dt:Num_period*T1-dt;
x = zeros(length(t), size(x0,1)); 
for j = 1:size(x0,1)
    x(:,j) = interp1(tt, xx(:,j), t, 'spline'); 
end
figure(11)
plot(x(:,1),x(:,2));

X = fft(x)/length(x(:,1)); X(2:end,:)=X(2:end,:)*2;
na = length(x(:,1)); ff = (0:na-1)'/na*2*pi/dt;
[pks,locs] = findpeaks(log(abs(X(1:end/2,1))));
[pks_sorted, idx_sorted] = sort(pks, 'descend');
keep_idx = pks_sorted >= (pks_sorted(1) - 12);  
pks_final = pks_sorted(keep_idx);
loc_final = locs(idx_sorted(keep_idx));
loc_final = [1;loc_final];
ff_final = ff(loc_final,:)/(OM0(1));
X_final  = X(loc_final,:);
X_MF1 = solve_HOrder(X_final,ff_final,rho,5,128);
X_MF2 = X_MF1(:,1:end-1);
%----------------------------------------------------------
figure(12)
plot(ff/(OM0(1)),log(abs(X(:,1))))
xlabel('normalized f')
ylabel('log(abs(X))')
xlim([0,4])
k_h = X_MF1(:,1)';
Q   = X_MF1(:,2:end);
theta = (OM0*t)';
H_IDF = exp(1i*theta*k_h);
xx = real(H_IDF*Q);
figure(13)
plot(x(:,1),x(:,2))
hold on
plot(xx(:,1),xx(:,2),'--')
xlabel('x')
ylabel('dx')
legend('TI','TI-FFT')
%----------------------------------------------------------
%%
prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',1, 'e',1, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',false,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',2^5,'L',2^8,'N_Ly',1000,'Omi',1);
prob1 = coco_set(prob1,'cont','PtMX',1000);
prob1 = coco_set(prob1,'cont','NSV',10);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.1);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB'} }, ...   
    'k',       { {11} }, ...        
    'P',       { {128} }, ...    
    'S',       { {2^5} }, ...
    'S_wave',  { {2^8} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname, p0);
bd1 = coco(prob1, 'p64n', [], 0, {'k_OM1' 'r' 'b' 'f'}, [2 5]);

% bd1 = coco_bd_read('qp1');
max_x1    = coco_bd_col(bd1, 'max.x');
min_x1    = coco_bd_col(bd1, 'min.x');
Om_mode1  = coco_bd_col(bd1, 'Om');
Stab      = coco_bd_col(bd1, 'Stab');
x_ht      = coco_bd_col(bd1, 'HT_x');
u_ht      = coco_bd_col(bd1, 'HT_u');
A = find(~Stab);

% figure(1)
% hold on
% for ii = 0:127
% plot(Om_mode1, x_ht(2^8*ii+1,:),'m-');
% % plot(Om_mode1(A), x_ht(2^8*ii+1,A), 'yo');
% end
% 
% xlabel('\omega_1');
% ylabel('x');
% grid on;



figure(2)
plot(x(:,1),x(:,2),'r-');
hold on
plot(x_ht(:,end),u_ht(:,end),'b--');




