clearvars;
close all;
clc;
% system parameters
OM0 = 1;
pname = {'Q'};
p0 = 6;
%%
% initial condition
% ---- time intergration method
x0 = zeros(6,1)+0.1;
Num_period = 2^12;
T1 = 2*pi/(OM0(1));
S_period = 2^12;
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
X = fft(x)/length(x(:,1)); X(2:end,:)=X(2:end,:)*2;
na = length(x(:,1)); ff = (0:na-1)'/na*2*pi/dt;
[pks,locs] = findpeaks(log(abs(X(1:end/2,1))));
[pks_sorted, idx_sorted] = sort(pks, 'descend');
keep_idx = pks_sorted >= (pks_sorted(1) - 8);  
pks_final = pks_sorted(keep_idx);
loc_final = locs(idx_sorted(keep_idx));
ff_final = ff(loc_final,:)/(OM0(1));
X_final  = X(loc_final,:);

figure(1)
plot(ff/(OM0(1)),log(abs(X(:,1))))
xlabel('normalized f')
ylabel('log(abs(X))')
xlim([0,5])

rho = [ff_final(1)];
X_MF1 = solve_HOrder(X_final,ff_final,rho,[5],[1]);
X_MF2 = X_MF1(:,1:4);
% ---- analytical method
OM0 = rho;
k_h = X_MF1(:,1)';
Q   = X_MF1(:,2:end);
theta = (OM0*t)';
H_IDF = exp(1i*theta*k_h);
xx = real(H_IDF*Q);

figure(2)
plot(x(:,2),x(:,5))
hold on
plot(xx(:,2),xx(:,5),'--')
xlabel('q2')
ylabel('dq2')
legend('TI','TI-FFT')

prob = coco_prob();
ode_opts = struct('n',3, 'd',1, 'e',0, 'isfwithu',false);
prob = coco_set(prob,'cont','PtMX',1000);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',2,'Nperiod',10,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',500,'Omi',1);
prob = coco_set(prob,'cont','h_min',1e-3,'h_max',0.1);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'           , ...     
    'type',    { {'HB'} }, ...   
    'k',       { {1:1:9} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^5} }, ...
    'S_wave',  { {2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Airfoil1_F, @Airfoil1_DFDX, [], @Airfoil1_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF2, OM0, pname, p0);

bd = coco(prob, 'qp', [], 1, { 'Q' 'u_om1'}, [4.75 10]);
% close(11111)

max_x    = coco_bd_col(bd, 'Max_x');
Om_mode  = coco_bd_col(bd, 'Om');
Stab     = coco_bd_col(bd, 'Stab');
Ly       = coco_bd_col(bd, 'Ly');
p        = coco_bd_col(bd, 'p');

HT_x = coco_bd_col(bd,'HT_x');
HT_u = coco_bd_col(bd,'HT_u');
T_x = coco_bd_col(bd,'T_x');
T_u = coco_bd_col(bd,'T_u');

MF_x = coco_bd_col(bd,'MF_x');
MF_f = coco_bd_col(bd,'MF_f');

figure(6)
plot(ff/(OM0(1)),log(abs(X(:,1))),'k-')
hold on
plot(MF_f{1}(:,1),log(abs(2*MF_x{1}(:,1))),'ro')
xlabel('normalized f')
ylabel('log(abs(X))')
xlim([0,6])

figure(7)
plot(x(:,2),x(:,5))
hold on
plot(HT_x(:,1,1),HT_u(:,1,1),'.')
xlabel('q2')
ylabel('dq2')
