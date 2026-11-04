clearvars;
close all;
clc;

pname = {'Q'};
p0 = 7.11609946946619; 
OM0 = [0.554933536226677
0.138636696133489];

% load('IP_QP.mat');
% X_MF = [IP_QP.location{1}',transpose(IP_QP.amplitude{1})];
% for ii = 2:numel(IP_QP.location)
% X_MF = [X_MF;IP_QP.location{ii}',transpose(IP_QP.amplitude{ii})];
% end

X_MF = my_XMF;

prob = coco_prob();
ode_opts = struct('n',3, 'd',2, 'e',0, 'isfwithu',false);
prob = coco_set(prob,'cont','PtMX',1000);
amp_opts = struct('isamp',false,'isuamp',true,'outdof',2,'Nperiod',10,'Nsampling',2^8);
sta_opts = struct('issta',false,'S',2^6,'L',2^8,'N_Ly',500,'Omt',[1,2]);
par_opts    = struct('parallel',true,'ncores',6);
prob = coco_set(prob,'cont','NSV',1);
prob = coco_set(prob,'cont','h_min',1e-2,'h_max',1e-1);
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts,'par_opts',par_opts); 
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @Newmark, ...
    'isvar'    ,  false, ...
    'Omt',        1, ...
    'M_SHT',      2, ...
    'Nsteps',     2^10,...
    'k',          { {11,5} }, ...        
    'P',          { {1,1} }, ...    
    'S',          { {2^5,2^5} }, ...
    'S_wave',     { {2^8,2^8} });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Airfoil1_F, @Airfoil1_DFDX, [], @Airfoil1_L}; 
prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF, OM0, pname, p0);
bd = coco(prob, 'qp', [], 1, { 'Q' 'u_om1' 'u_om2'}, [p0 9]);

% bd = coco_bd_read('qp');
Om1  = coco_bd_col(bd, 'u_om1');
Om2  = coco_bd_col(bd, 'u_om2');
figure(2)
plot(Om1,Om2);


Om = [Om1;Om2];