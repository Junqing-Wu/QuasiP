clearvars;
close all;
clc;

num_cont = 100;

e1      = 2;
e2      = 1;
e3      = 0.5;

Fex1.location{1} = [1;0;0];
Fex1.amplitude{1} = e1;
Fex1.location{2} = [0;1;0];
Fex1.amplitude{2} = e2;
Fex1.location{3} = [0;0;1];
Fex1.amplitude{3} = e3;

om1 = 10;
rho = [1;1/5^0.5;1/11^0.5];
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

%% ===================== 任务：SHT =====================

prob = coco_prob();
ode_opts    = struct('n',1, 'd',3, 'e',3, 'isfwithu',true);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',1,'Nperiod',1,'Nsampling',2^8);
par_opts    = struct('parallel',true,'ncores',12);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'par_opts',par_opts); 
prob = coco_set(prob,'cont','PtMX',num_cont);
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @ode45, ...
    'isvar',      true ,...
    'Omt',        1, ...
    'M_SHT',      1, ...
    'Nsteps',     2^8, ...   
    'k',       { {5,5,5} }, ...        
    'P',       { {1,1,1} }, ...    
    'S',       { {2^4,2^4,2^4} }, ...
    'S_wave',  { {2^6,2^6,2^6} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@VDPDwithe3_F, @VDPD_DFDX, [], @VDPD_L, @VDPDwithe3_DFDTH}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_3e};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');

bd3 = coco(prob, 'sht3', [], 1, {'k_OM1' 'k_OM2' 'k_OM3'}, [0.2 om1]);
close(11111)

max_x1    = coco_bd_col(bd3, 'Max_x');
Om_mode  = coco_bd_col(bd3, 'Om');

figure(1)
hold on
plot(Om_mode(:,1), max_x1,'-','DisplayName','HB+HB+HB');
xlabel('\omega_1');
ylabel('max(x)');
grid on;


function [data, y, J] = fc_3e(prob, data, u)
rho = [1;1/5^0.5;1/11^0.5];
Om1 = u(1);
Om2 = u(2);
Om3 = u(3);

y = [rho(2)* Om1 - Om2;
     rho(3)* Om1 - Om3];
J = [rho(2), -1,  0;
     rho(3),  0, -1];
end



