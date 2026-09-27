clearvars;
close all;
clc;
% system parameters
rho = 1;
Om1 = 0.9;
OM0 = rho*Om1;
pname = {'BB' 'knl' 'AA'};
p0 = [0.4;3;0.3];
%----------------------------------------------------------
% initial condition
% ---- time intergration method
x0 = zeros(4,1);
Num_period = 2^10;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = false;[t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@NES_F,@NES_L,isfwithu1), tspan, x0, opts);
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
X_MF1 = solve_HOrder(X_final,ff_final,rho,5,1);
X_MF2 = X_MF1(:,1:end-2);
%% ====================== 任务 1 ======================
prob1 = coco_prob();
ode_opts    = struct('n',2, 'd',1, 'e',1, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'cont','h_min',0.001,'h_max',0.01);

amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',1,'L',2^8,'N_Ly',[],'Omi',1);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'         , ...  
    'isvar'    ,true          , ...
    'method',  { {'HB'} }, ...  
    'k',       { {5} }      , ...        
    'P',       { {1} }      , ...    
    'S',       { {2^6} }  , ...
    'S_wave',  { {2^8} }  , ...
    'OmOd' ,   1);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@NES_F, @NES_DFDX, @NES_DFDP, @NES_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
bd1 = coco(prob1, 'qp1', [], 1, {'k_OM1' 'BB' 'knl' 'AA'}, [0.9 1.1]);
% close(11111)

% bd1 = coco_bd_read('qp1');
max_x1    = coco_bd_col(bd1, 'Max_x');
Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1);

figure(1)
plot(Om_mode1, max_x1(:,1),'k--');
hold on
plot(Om_mode1(A11), max_x1(A11,1), 'ko');
xlabel('\omega_1');
ylabel('max(x)');
grid on;

%% continuation of QP from TR bifurcation
lab  = coco_bd_labs(bd1, 'TR');
lab  = max(lab);
prob2 = coco_prob();
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
stable_opts = struct('issta',true,'S',2^5,'L',2^8,'N_Ly',500,'Omi',[1,2]);
prob2 = coco_set(prob2,'QuasiP','amp_opts',amp_opts,'stable_opts',stable_opts); 

% PHI_INITIAL = struct( 'parametrization_tori', 'MT', ...
%     'type',       'SHT'           , ...     
%     'method',     'MF' , ... % HT or MF
%     'ODEsolver',  @ode45, ...
%     'isvar'    ,  true, ...
%     'Omt',        1, ...
%     'M_SHT',      1, ...
%     'Nsteps',     2^8,...
%     'k',          { {5,5} }, ...        
%     'P',          { {1,1} }, ...    
%     'S',          { {2^5,2^5} }, ...
%     'S_wave',     { {2^8,2^8} });
% prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
% prob2 = ode_TR2QuasiP(prob2, '', 'qp1',  lab , @NES_DFDTH);
% bd2 = coco(prob2, 'qp2', [], 1, { 'k_OM1' 'u_om2' 'BB' 'knl' 'AA'}, [0.971 1.035]);

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'         , ...  
    'isvar'    ,true          , ...
    'method',  { {'HB','HB'} }, ...  
    'k',       { {5,5} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^5,2^5} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
prob2 = ode_TR2QuasiP(prob2, '', 'qp1',  lab );
bd2 = coco(prob2, 'qp2', [], 1, { 'k_OM1' 'u_om2' 'BB' 'knl' 'AA'}, [0.971 1.035]);



