clearvars;
close all;
clc;

rho = [1;1/2^0.5];
Om1 = 1;
OM0 = rho*Om1;
X_MF2 = [1,0,1;
         0,1,1;
         1,1,1];

prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',2, 'e',0, 'isfwithu',false);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts);
prob1 = coco_set(prob1,'cont','PtMX',3000);
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'         , ...  
    'isvar'    ,true          , ...
    'type',    { {'HB','HB'} }, ...  
    'k',       { {3,3} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^6,2^6} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@Duffing_F, @Duffing_DFDX};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0);