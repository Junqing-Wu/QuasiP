clearvars;
close all;
clc;


num_cont = 2;
om1 = 2.65;
rho = [1;1/5^0.5;1/11^0.5;1/33^0.5];
OM0  = om1*rho;
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
isfwithu1 = true;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@VDPDwithe4_F,@VDPD_L,isfwithu1), tspan, x0, opts);

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
X_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5,5,5],[1,1,1,1]);
X_MF2 = X_MF1(:,1:end-1);
% ----------------------------------------------------------

prob = coco_prob();
ode_opts = struct('n',1, 'd',4, 'e',4, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',10000,'Nsampling',2^2);
prob = coco_set(prob,'cont','PtMX',num_cont);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB','HB','HB'} }, ...   
    'k',       { {5,3,3,3} }, ...        
    'P',       { {1,1,1,1} }, ...    
    'S',       { {2^5,2^5,2^5,2^5} }, ...
    'S_wave',  { {2^6,2^6,2^6,2^6} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe4_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF2, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_4e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');

bd4_hbhbhbhb = coco(prob, 'qp4', [], 0, {'k_OM1' 'k_OM2' 'k_OM3' 'k_OM4'}, [0.2 om1]);

max_x1_hbhbhbhb    = coco_bd_col(bd4_hbhbhbhb, 'Max_x');
Om_mode_hbhbhbhb  = coco_bd_col(bd4_hbhbhbhb, 'Om');
HT_x  = coco_bd_col(bd4_hbhbhbhb, 'HT_x');
HT_u  = coco_bd_col(bd4_hbhbhbhb, 'HT_u');
MF_x  = coco_bd_col(bd4_hbhbhbhb, 'MF_x');
MF_f  = coco_bd_col(bd4_hbhbhbhb, 'MF_f');
T_x  = coco_bd_col(bd4_hbhbhbhb, 'T_x');
T_u  = coco_bd_col(bd4_hbhbhbhb, 'T_u');

figure(1)
plot(ff/(OM0(1)),log(abs(X(:,1))),'k-')
hold on
plot(MF_f(:,1),log(abs(2*MF_x)),'ro')
xlabel('normalized f')
ylabel('log(abs(X))')
title('$x(\tau_1,\tau_2,\tau_3,\tau_4)$','Interpreter','latex');
xlim([0,6])

HT_x_4 = reshape(HT_x,2^6,2^6,2^6,2^6);
HT_u_4 = reshape(HT_u,2^6,2^6,2^6,2^6);

figure(2)
%% 221: Full 4-torus solution
subplot(221)
plot(HT_x,HT_u,'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_4(1,:,:,:),[],1),reshape(HT_u_4(1,:,:,:),[],1),'b-')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(\tau_1,\tau_2,\tau_3,\tau_4)$','Interpreter','latex');
grid on;

%% 222: Poincare section 3D → 2D
subplot(222)
plot(reshape(HT_x_4(1,:,:,:),[],1),reshape(HT_u_4(1,:,:,:),[],1),'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_4(1,1,:,:),[],1),reshape(HT_u_4(1,1,:,:),[],1),'b-')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(0,\tau_2,\tau_3,\tau_4)$','Interpreter','latex');
grid on;

%% 223: Poincare section 2D → 1D
subplot(223)
plot(reshape(HT_x_4(1,1,:,:),[],1),reshape(HT_u_4(1,1,:,:),[],1),'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_4(1,1,1,:),[],1),reshape(HT_u_4(1,1,1,:),[],1),'b-')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(0,0,\tau_3,\tau_4)$','Interpreter','latex');
grid on;

%% 224: Poincare section 1D → periodic
subplot(224)
plot(reshape(HT_x_4(1,1,1,:),[],1),reshape(HT_u_4(1,1,1,:),[],1),'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_4(1,1,1,1),[],1),reshape(HT_u_4(1,1,1,1),[],1),'bo')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(0,0,0,\tau_4)$','Interpreter','latex');
grid on;

%%

num_cont = 2;
om1 = 2.65;
rho = [1;1/5^0.5;1/11^0.5]; 
OM0  = om1*rho;

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
isfwithu1 = true;
[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@VDPDwithe3_F,@VDPD_L,isfwithu1), tspan, x0, opts);

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
X_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5,5],[1,1,1]);  
X_MF2 = X_MF1(:,1:end-1);
% ----------------------------------------------------------

prob = coco_prob();
ode_opts = struct('n',1, 'd',3, 'e',3, 'isfwithu',true); 
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',10000,'Nsampling',2^2);
prob = coco_set(prob,'cont','PtMX',num_cont);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB','HB'} }, ...  
    'k',       { {5,5,5} }, ...           
    'P',       { {1,1,1} }, ...          
    'S',       { {2^5,2^5,2^5} }, ...     
    'S_wave',  { {2^6,2^6,2^6} });        
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@VDPDwithe3_F, @VDPD_DFDX, [], @VDPD_L};  
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF2, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_3e};  
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');

bd3_hbhbhb = coco(prob, 'qp3', [], 0, {'k_OM1' 'k_OM2' 'k_OM3'}, [0.2 om1]);

max_x1_hbhbhb    = coco_bd_col(bd3_hbhbhb, 'Max_x');
Om_mode_hbhbhb  = coco_bd_col(bd3_hbhbhb, 'Om');
HT_x  = coco_bd_col(bd3_hbhbhb, 'HT_x');
HT_u  = coco_bd_col(bd3_hbhbhb, 'HT_u');
MF_x  = coco_bd_col(bd3_hbhbhb, 'MF_x');
MF_f  = coco_bd_col(bd3_hbhbhb, 'MF_f');
T_x  = coco_bd_col(bd3_hbhbhb, 'T_x');
T_u  = coco_bd_col(bd3_hbhbhb, 'T_u');

figure(11)
plot(ff/(OM0(1)),log(abs(X(:,1))),'k-')
hold on
plot(MF_f(:,1),log(abs(2*MF_x)),'ro')
xlabel('normalized f')
ylabel('log(abs(X))')
title('$x(\tau_1,\tau_2,\tau_3)$','Interpreter','latex');
xlim([0,6])

HT_x_3 = reshape(HT_x,2^6,2^6,2^6);
HT_u_3 = reshape(HT_u,2^6,2^6,2^6);

figure(22)
%% 221: Full 3-torus solution
subplot(221)
plot(HT_x,HT_u,'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_3(1,:,:),[],1),reshape(HT_u_3(1,:,:),[],1),'b-')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(\tau_1,\tau_2,\tau_3)$','Interpreter','latex');
grid on;

%% 222: Poincare section 3D → 2D
subplot(222)
plot(reshape(HT_x_3(1,:,:),[],1),reshape(HT_u_3(1,:,:),[],1),'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_3(1,1,:),[],1),reshape(HT_u_3(1,1,:),[],1),'b-')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(0,\tau_2,\tau_3)$','Interpreter','latex');
grid on;

%% 223: Poincare section 2D → 1D
subplot(223)
plot(reshape(HT_x_3(1,1,:),[],1),reshape(HT_u_3(1,1,:),[],1),'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_3(1,1,1),[],1),reshape(HT_u_3(1,1,1),[],1),'bo')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(0,0,\tau_3)$','Interpreter','latex');
grid on;

subplot(224)
axis off

num_cont = 2;
om1 = 2.65;
rho = [1;1/5^0.5]; 
OM0  = om1*rho;

x0 = zeros(2,1);
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = true;
[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@VDPDwithe2_F,@VDPD_L,isfwithu1), tspan, x0, opts);

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
X_MF1 = solve_HOrder(X_final,ff_final,rho,[5,5],[1,1]);  
X_MF2 = X_MF1(:,1:end-1);

prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true); 
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',10000,'Nsampling',2^2);
prob = coco_set(prob,'cont','PtMX',num_cont);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB'} }, ...  
    'k',       { {5,5} }, ...           
    'P',       { {1,1} }, ...          
    'S',       { {2^5,2^5} }, ...     
    'S_wave',  { {2^6,2^6} });        
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L};  
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF2, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};  
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');

bd2_hbhb = coco(prob, 'qp2', [], 0, {'k_OM1' 'k_OM2'}, [0.2 om1]);

max_x1_hbhb    = coco_bd_col(bd2_hbhb, 'Max_x');
Om_mode_hbhb  = coco_bd_col(bd2_hbhb, 'Om');
HT_x  = coco_bd_col(bd2_hbhb, 'HT_x');
HT_u  = coco_bd_col(bd2_hbhb, 'HT_u');
MF_x  = coco_bd_col(bd2_hbhb, 'MF_x');
MF_f  = coco_bd_col(bd2_hbhb, 'MF_f');
T_x  = coco_bd_col(bd2_hbhb, 'T_x');
T_u  = coco_bd_col(bd2_hbhb, 'T_u');

figure(111)
plot(ff/(OM0(1)),log(abs(X(:,1))),'k-')
hold on
plot(MF_f(:,1),log(abs(2*MF_x)),'ro')
xlabel('normalized f')
ylabel('log(abs(X))')
title('$x(\tau_1,\tau_2)$','Interpreter','latex');
xlim([0,6])

HT_x_2 = reshape(HT_x,2^6,2^6);
HT_u_2 = reshape(HT_u,2^6,2^6);

figure(222)
subplot(221)
plot(HT_x,HT_u,'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_2(1,:),[],1),reshape(HT_u_2(1,:),[],1),'b-')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(\tau_1,\tau_2)$','Interpreter','latex');
grid on;

subplot(222)
plot(reshape(HT_x_2(1,:),[],1),reshape(HT_u_2(1,:),[],1),'--','Color',[0.5 0.5 0.5])
hold on
plot(reshape(HT_x_2(1,1),[],1),reshape(HT_u_2(1,1),[],1),'bo')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(0,\tau_2)$','Interpreter','latex');
grid on;

subplot(223)
axis off
subplot(224)
axis off


num_cont = 2;
om1 = 2.65;
rho = 1; 
OM0  = om1*rho;

x0 = zeros(2,1);
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = true;
[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,[],OM0,@VDPDwithe1_F,@VDPD_L,isfwithu1), tspan, x0, opts);

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
X_MF1 = solve_HOrder(X_final,ff_final,rho,5,1);  
X_MF2 = X_MF1(:,1:end-1);

prob = coco_prob();
ode_opts = struct('n',1, 'd',1, 'e',1, 'isfwithu',true); 
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',10000,'Nsampling',2^2);
prob = coco_set(prob,'cont','PtMX',num_cont);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB'} }, ...  
    'k',       { {5} }, ...           
    'P',       { {1} }, ...          
    'S',       { {2^5} }, ...     
    'S_wave',  { {2^6} });        
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@VDPDwithe1_F, @VDPD_DFDX, [], @VDPD_L};  
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF2, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');

bd1_hb = coco(prob, 'qp1', [], 0, {'k_OM1'}, [0.2 om1]);

max_x1_hb    = coco_bd_col(bd1_hb, 'Max_x');
Om_mode_hb  = coco_bd_col(bd1_hb, 'Om');
HT_x  = coco_bd_col(bd1_hb, 'HT_x');
HT_u  = coco_bd_col(bd1_hb, 'HT_u');
MF_x  = coco_bd_col(bd1_hb, 'MF_x');
MF_f  = coco_bd_col(bd1_hb, 'MF_f');
T_x  = coco_bd_col(bd1_hb, 'T_x');
T_u  = coco_bd_col(bd1_hb, 'T_u');

figure(1111)
plot(ff/(OM0(1)),log(abs(X(:,1))),'k-')
hold on
plot(MF_f(:,1),log(abs(2*MF_x)),'ro')
xlabel('normalized f')
ylabel('log(abs(X))')
title('$x(\tau_1)$','Interpreter','latex');
xlim([0,6])

figure(2222)
subplot(221)
plot(HT_x,HT_u,'--','Color',[0.5 0.5 0.5])
hold on
plot(HT_x(1),HT_u(1),'bo')
xlabel('$x$','Interpreter','latex'); 
ylabel('$\dot{x}$','Interpreter','latex');
title('$x(\tau_1)$','Interpreter','latex');
grid on;

subplot(222)
axis off
subplot(223)
axis off
subplot(224)
axis off


function [data, y, J] = fc_4e(prob, data, u)

rho = [1;1/5^0.5;1/11^0.5;1/33^0.5];
Om1 = u(1);
Om2 = u(2);
Om3 = u(3);
Om4 = u(4);

y = [rho(2)* Om1 - Om2;
     rho(3)* Om1 - Om3;
     rho(4)* Om1 - Om4];
J = [rho(2), -1,  0,  0;
     rho(3),  0, -1,  0;
     rho(4),  0,  0, -1];
end

function [data, y, J] = fc_3e(prob, data, u)

rho = [1;1/5^0.5;1/11^0.5];
Om1 = u(1);
Om2 = u(2);
Om3 = u(3);

y = [rho(2)* Om1 - Om2;
     rho(3)* Om1 - Om3];
J = [rho(2), -1,  0;
     rho(3),  0, -1,];
end

function [data, y, J] = fc_2e(prob, data, u)

rho = [1;1/5^0.5];
Om1 = u(1);
Om2 = u(2);

y = [rho(2)* Om1 - Om2];
J = [rho(2), -1];
end