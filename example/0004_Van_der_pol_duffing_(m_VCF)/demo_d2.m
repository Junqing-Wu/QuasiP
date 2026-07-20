clearvars;
close all;
clc;


num_cont = 1000;

e1      = 2;
e2      = 1;
Fex1.location{1} = [1;0];
Fex1.amplitude{1} = e1;
Fex1.location{2} = [0;1];
Fex1.amplitude{2} = e2;

om1 = 9;
rho = [1;1/5^0.5];
OM0  = om1*rho;

[~,~,A,B,C] = VDPD_L([]);

p = length(OM0);
for ii = 1:p
    IP.amplitude{ii} = ( -(OM0(ii))^2*A + 1i*OM0(ii)*B + C )\ ...
       Fex1.amplitude{ii} ;
    IP.location{ii} = Fex1.location{ii};
end

X_MF = [IP.location{1}',transpose(IP.amplitude{1})];
for ii = 2:numel(IP.location)
X_MF = [X_MF;IP.location{ii}',transpose(IP.amplitude{ii})];
end

%% ===================== 任务1：HB + HB =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB'} }, ...   
    'k',       { {5,5} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} }, ...
    'OmOd',    [2,1]);
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_hbhb = coco(prob, 'qphb+hb', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_hbhb    = coco_bd_col(bd2_hbhb, 'Max_x');
Om_mode_hbhb  = coco_bd_col(bd2_hbhb, 'Om');

figure(1)
hold on
plot(Om_mode_hbhb(:,1), max_x1_hbhb,'-','DisplayName','HB+HB');
xlabel('\omega_1');
ylabel('max(x)');
grid on;

%% ===================== 任务2：HB + CO =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','CO'} }, ...   
    'k',       { {5,4} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_hbco = coco(prob, 'qphb+co', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_hbco    = coco_bd_col(bd2_hbco, 'Max_x');
Om_mode_hbco  = coco_bd_col(bd2_hbco, 'Om');

figure(1)
hold on
plot(Om_mode_hbco(:,1), max_x1_hbco,'-','DisplayName','HB+CO');

%% ===================== 任务3：HB + FD =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','FD'} }, ...   
    'k',       { {5,-2:1:1} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_hbfd = coco(prob, 'qphb+fd', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_hbfd    = coco_bd_col(bd2_hbfd, 'Max_x');
Om_mode_hbfd  = coco_bd_col(bd2_hbfd, 'Om');

figure(1)
hold on
plot(Om_mode_hbfd(:,1), max_x1_hbfd,'-','DisplayName','HB+FD');

%% ===================== 任务4：CO + HB =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'CO','HB'} }, ...   
    'k',       { {4,5} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_cohb = coco(prob, 'qpco+hb', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_cohb    = coco_bd_col(bd2_cohb, 'Max_x');
Om_mode_cohb  = coco_bd_col(bd2_cohb, 'Om');

figure(1)
hold on
plot(Om_mode_cohb(:,1), max_x1_cohb,'-','DisplayName','CO+HB');

%% ===================== 任务5：CO + CO =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'CO','CO'} }, ...   
    'k',       { {4,4} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_coco = coco(prob, 'qpco+co', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_coco    = coco_bd_col(bd2_coco, 'Max_x');
Om_mode_coco  = coco_bd_col(bd2_coco, 'Om');

figure(1)
hold on
plot(Om_mode_coco(:,1), max_x1_coco,'-','DisplayName','CO+CO'); % 修复原始笔误

%% ===================== 任务6：CO + FD =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'CO','FD'} }, ...   
    'k',       { {4,-2:1:1} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_cofd = coco(prob, 'qpco+fd', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_cofd    = coco_bd_col(bd2_cofd, 'Max_x');
Om_mode_cofd  = coco_bd_col(bd2_cofd, 'Om');

figure(1)
hold on
plot(Om_mode_cofd(:,1), max_x1_cofd,'-','DisplayName','CO+FD');

%% ===================== 任务7：FD + HB =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'FD','HB'} }, ...   
    'k',       { {-2:1:1,5} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_fdhb = coco(prob, 'qpfd+hb', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_fdhb    = coco_bd_col(bd2_fdhb, 'Max_x');
Om_mode_fdhb  = coco_bd_col(bd2_fdhb, 'Om');

figure(1)
hold on
plot(Om_mode_fdhb(:,1), max_x1_fdhb,'-','DisplayName','FD+HB');

%% ===================== 任务8：FD + CO =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'FD','CO'} }, ...   
    'k',       { {-2:1:1,4} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_fdco = coco(prob, 'qpfd+co', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_fdco    = coco_bd_col(bd2_fdco, 'Max_x');
Om_mode_fdco  = coco_bd_col(bd2_fdco, 'Om');

figure(1)
hold on
plot(Om_mode_fdco(:,1), max_x1_fdco,'-','DisplayName','FD+CO');

%% ===================== 任务9：FD + FD =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'FD','FD'} }, ...   
    'k',       { {-2:1:1,-2:1:1} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe2_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd2_fdfd = coco(prob, 'qpfd+fd', [], 1, {'k_OM1' 'k_OM2'}, [0.2 om1]);
close(11111)

max_x1_fdfd    = coco_bd_col(bd2_fdfd, 'Max_x');
Om_mode_fdfd  = coco_bd_col(bd2_fdfd, 'Om');

figure(1)
hold on
plot(Om_mode_fdfd(:,1), max_x1_fdfd,'-','DisplayName','FD+FD');
legend('Location','best');
grid on;



function [data, y, J] = fc_2e(prob, data, u)
rho = [1;1/5^0.5];
Om1 = u(1);
Om2 = u(2);

y = rho(2)* Om1 - Om2;
J = [rho(2), -1];
end