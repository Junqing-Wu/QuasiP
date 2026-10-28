%% ====================== tor44 ======================
% system parameters
rho = [1;1/2^0.5];
Om1 = 2.71;%2.71----2.75----2.85
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% initial condition
% ---- time intergration method
bd_tor44   = coco_bd_read('qp515_tor44');
MF_f       = coco_bd_col(bd_tor44,'MF_f');
MF_x       = coco_bd_col(bd_tor44,'MF_x');
X_MF2_tor44 = [MF_f(:,3:4),MF_x];
%----------------------------------------------------------
%% ====================== tor44 COCO ======================
prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.5);
stable_opts = struct('issta',true,'S',2^5,'L',2^8,'N_Ly',2000,'Omi',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {4,4} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^7,2^7} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob1, 'qp514_tor44_a', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [2.3236 2.71]);


prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.25);
stable_opts = struct('issta',true,'S',2^5,'L',2^8,'N_Ly',2000,'Omi',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {4,4} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^7,2^7} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob1, 'qp514_tor44_b', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [2.71 2.72]);


function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end
