clc;
clear;

bd2 = coco_bd_read('qp2');

SMF_x = coco_bd_col(bd2,'SX.MF_x');
SMF_f = coco_bd_col(bd2,'SX.MF_f');
SX_MF  = [SMF_f{1}(:,3:4),SMF_x{1}];

USMF_x = coco_bd_col(bd2,'USX.MF_x');
USMF_f = coco_bd_col(bd2,'USX.MF_f');
USX_MF  = [USMF_f{1}(:,3:4),USMF_x{1}];

SOM = coco_bd_col(bd2,'SX.Om');
USOM = coco_bd_col(bd2,'USX.Om');

SOM0    = SOM(1,:)';
USOM0   = USOM(1,:)';

pname = {'mu','OM'};
p0 = [0.2;SOM0(1)];

prob3 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',1, 'isfwithu',true);
amp_opts = struct('isamp',false,'isuamp',true,'outdof',[1 2],'Nperiod',50,'Nsampling',2^8);
par_opts = struct('parallel',true,'ncores',12);

prob3 = coco_set(prob3,'cont','PtMX',3000);
prob3 = coco_set(prob3,'cont','h_min',1e-3,'h_max',0.5);
prob3 = coco_set(prob3,'cont','NSV',1);
prob3 = coco_set(prob3,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'par_opts',par_opts); 
% prob3 = coco_set(prob3, 'cont','PtMX',25);
% prob3 = coco_set(prob3, 'corr','ItMX',50);
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @ode45, ...
    'isvar'    ,  true, ...
    'Omt',        2, ...
    'M_SHT',      1, ...
    'Nsteps',     2^9,... 
    'k',       { {9,9} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^8,2^8} }, ...
    'S_wave',  { {2^8,2^8} });
prob3 = coco_set(prob3, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@RS1_F, @RS1_DFDX, @RS1_DFDP, @RS1_L, @RS1_DFDTH};
prob3 = ode_isol2QuasiP(prob3, 'SX', funcs{:}, SX_MF, SOM0, pname , p0);
[SXdata, SXuidx] = coco_get_func_data(prob3, 'SX.QuasiP', 'data', 'uidx');
prob3 = coco_add_glue(prob3,'SX.OMp2',SXuidx(SXdata.kOM_idx),SXuidx(SXdata.p_idx(2)));

prob3 = ode_isol2QuasiP(prob3, 'USX', funcs{:}, USX_MF, USOM0, pname , p0);
[USXdata, USXuidx] = coco_get_func_data(prob3, 'USX.QuasiP', 'data', 'uidx');
prob3 = coco_add_glue(prob3,'USX.OMp2',USXuidx(USXdata.kOM_idx),USXuidx(USXdata.p_idx(2)));

prob3 = coco_add_glue(prob3,'pname1',SXuidx(SXdata.p_idx(1)),USXuidx(USXdata.p_idx(1)));
prob3 = coco_add_glue(prob3,'pname2',SXuidx(SXdata.p_idx(2)),USXuidx(USXdata.p_idx(2)));
% prob3 = coco_add_glue(prob3,'KOM1',  SXuidx(SXdata.kOM_idx),USXuidx(USXdata.kOM_idx));

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

% plot(SOM1,log(Smu-0.1),'r--')


function [data y J] = Dom(prob,data,u)
SOm2 = u(1);
USOm2 = u(2);
y = SOm2 - USOm2;
J = [1 -1];
end

