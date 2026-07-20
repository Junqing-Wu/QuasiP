clearvars;
close all;
clc;

% 系统参数
e1      = 2;
Fex1.location{1} = 1;
Fex1.amplitude{1} = e1;

om1 = 10;
rho = 1;
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

%% ===================== 任务1：HB =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',1, 'e',1, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',1000);

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB'} }, ...   
    'k',       { {3} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^5} }, ...
    'S_wave',  { {2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@VDPDwithe1_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
bd1_hb = coco(prob, 'qp1_hb', [], 1, {'k_OM1'}, [0.2 om1]);
close(11111)

max_x1_hb    = coco_bd_col(bd1_hb, 'Max_x');
Om_mode_hb  = coco_bd_col(bd1_hb, 'Om');

figure(1)
plot(Om_mode_hb, max_x1_hb,'k--','DisplayName','HB');
xlabel('\omega_1');
ylabel('max(x)');
grid on;
hold on;

%% ===================== 任务2：CO =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',1, 'e',1, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',1000);

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'CO'} }, ...   
    'k',       { {4} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^5} }, ...
    'S_wave',  { {2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@VDPDwithe1_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
bd1_co = coco(prob, 'qp1_co', [], 1, {'k_OM1'}, [0.2 om1]);
close(11111)

max_x1_co    = coco_bd_col(bd1_co, 'Max_x');
Om_mode_co  = coco_bd_col(bd1_co, 'Om');
figure(1)
plot(Om_mode_co, max_x1_co,'r--','DisplayName','CO');

%% ===================== 任务3：FD =====================
prob = coco_prob();
ode_opts = struct('n',1, 'd',1, 'e',1, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
prob = coco_set(prob,'cont','PtMX',1000);

PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'FD'} }, ...   
    'k',       { {-2:1:1} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^5} }, ...
    'S_wave',  { {2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@VDPDwithe1_F, @VDPD_DFDX, [], @VDPD_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
bd1_fd = coco(prob, 'qp1_fd', [], 1, {'k_OM1'}, [0.2 om1]);
close(11111)

max_x1_fd    = coco_bd_col(bd1_fd, 'Max_x');
Om_mode_fd  = coco_bd_col(bd1_fd, 'Om');
figure(1)
plot(Om_mode_fd, max_x1_fd,'b.','DisplayName','FD');

legend('Location','best');
grid on;