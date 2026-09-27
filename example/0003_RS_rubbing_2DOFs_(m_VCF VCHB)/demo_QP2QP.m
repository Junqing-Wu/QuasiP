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
x0 = zeros(4,1)+0.0001;
Num_period = 2^9;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = true; [t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@RS1_F,@RS1_L,isfwithu1), tspan, x0, opts);

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
 
rho = ff_final(1:2);
X_MF = solve_HOrder(X_final,ff_final,rho,[5,5],[1,1]);
X_MF2 = X_MF(:,1:4);
OM0 = rho*Om1;

%% ==========================================

prob1 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',1, 'isfwithu',true);
amp_opts = struct('isamp',true,'outdof',[1 2],'Nperiod',1,'Nsampling',2^8);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB'} }, ...   
    'k',       { {7,7} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^8,2^8} }, ...
    'S_wave',  { {2^8,2^8} }, ...
    'OmOd',    [2,1]);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@RS1_F, @RS1_DFDX, @RS1_DFDP, @RS1_L};

prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname , p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
prob1 = coco_add_glue(prob1,'OMp2',data.kOM_idx,data.p_idx(2));
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

%% continuation of QP from previous solution
lab  = coco_bd_labs(bd1, 'EP');
lab1 = lab(1);
lab2 = lab(2);
prob2 = coco_prob();
prob2 = ode_QuasiP2QuasiP(prob2, 'SX', 'qp1', '', lab1);
[SXdata, SXuidx] = coco_get_func_data(prob2, 'SX.QuasiP', 'data', 'uidx');
prob2 = coco_add_glue(prob2,'SX.OMp2',SXuidx(SXdata.kOM_idx),SXuidx(SXdata.p_idx(2)));
%----------
prob2 = ode_QuasiP2QuasiP(prob2, 'USX', 'qp1', '', lab2);
[USXdata, USXuidx] = coco_get_func_data(prob2, 'USX.QuasiP', 'data', 'uidx');
prob2 = coco_add_glue(prob2,'USX.OMp2',USXuidx(USXdata.kOM_idx),USXuidx(USXdata.p_idx(2)));
%----------
prob2 = coco_add_glue(prob2,'KOM1',  SXuidx(SXdata.kOM_idx),USXuidx(USXdata.kOM_idx));
OM2_idx = [SXuidx(SXdata.uOM_idx),USXuidx(USXdata.uOM_idx)];
prob2 = coco_add_func(prob2,'Dom',@Dom,[],'inactive','Dom','uidx',OM2_idx,'F+dF');

Dom2  = om2(1)-om2(end);
bd2   = coco(prob2, 'qp2', [], 1, {'Dom' 'SX.u_om2' 'USX.u_om2' 'SX.k_OM1' 'USX.k_OM1' 'SX.OM' 'USX.OM' 'SX.mu' 'USX.mu' }, [0.002 Dom2]);
close(11111)

SOM1 = coco_bd_col(bd2,'SX.k_OM1');
Som2 = coco_bd_col(bd2,'SX.u_om2');
Smax_x = coco_bd_col(bd2, 'SX.Max_x');
USOM1 = coco_bd_col(bd2,'USX.k_OM1');
USom2 = coco_bd_col(bd2,'USX.u_om2');
USmax_x = coco_bd_col(bd2, 'USX.Max_x');

figure(1)
plot(OM1,om2,'r--');
hold on
plot(SOM1,Som2,'bo');
plot(USOM1,USom2,'k^');
xlabel('k_{OM1}');
ylabel('\omega_2');
grid on;

figure(2)
plot(OM1,max_x(:,1),'r--');
hold on
plot(SOM1,Smax_x(:,1),'bo');
plot(USOM1,USmax_x(:,1),'k^');
xlabel('k_{OM1}');
ylabel('max(x)');
grid on;


%% continuation of QP from previous solution
lab  = coco_bd_labs(bd2, 'EP');
lab  = max(lab);
prob3 = coco_prob();
prob3 = ode_QuasiP2QuasiP(prob3, 'SX', 'qp2', 'SX', lab);
[SXdata, SXuidx] = coco_get_func_data(prob3, 'SX.QuasiP', 'data', 'uidx');
prob3 = coco_add_glue(prob3,'SX.OMp2',SXuidx(SXdata.kOM_idx),SXuidx(SXdata.p_idx(2)));

prob3 = ode_QuasiP2QuasiP(prob3, 'USX', 'qp2', 'USX', lab);
[USXdata, USXuidx] = coco_get_func_data(prob3, 'USX.QuasiP', 'data', 'uidx');
prob3 = coco_add_glue(prob3,'USX.OMp2',USXuidx(USXdata.kOM_idx),USXuidx(USXdata.p_idx(2)));

prob3 = coco_add_glue(prob3,'pname1',SXuidx(SXdata.p_idx(1)),USXuidx(USXdata.p_idx(1)));
prob3 = coco_add_glue(prob3,'pname2',SXuidx(SXdata.p_idx(2)),USXuidx(USXdata.p_idx(2)));

OM2_idx = [SXuidx(SXdata.uOM_idx),USXuidx(USXdata.uOM_idx)];
prob3 = coco_add_func(prob3,'Dom',@Dom,[],'inactive','Dom','uidx',OM2_idx,'F+dF');

coco(prob3, 'qp3', [], 1, {'SX.mu' 'USX.mu' 'SX.u_om2' 'USX.u_om2' 'SX.k_OM1' 'USX.k_OM1' 'SX.OM' 'USX.OM' 'Dom'}, [0.1 0.4]);
close(11111)

bd3 = coco_bd_read('qp3');
SOM1 = coco_bd_col(bd3,'SX.k_OM1');
Smu = coco_bd_col(bd3,'SX.mu');

figure(1)
plot(SOM1,Smu,'r--')
xlim([0,1.2]);
ylim([0,0.4]);
xlabel('k_{OM1}');
ylabel('\mu');
grid on;


function [data y J] = Dom(prob,data,u)
SOm2 = u(1);
USOm2 = u(2);
y = SOm2 - USOm2;
J = [1 -1];
end
