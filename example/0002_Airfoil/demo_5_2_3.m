% DEMO_5_2_3  Construct an initial 2-torus problem from the upper TR point.
% Reference: Section 5.2, Eq. (35) initialization discussed with Fig. 14.
% Requires qp_d1 from demo_5_2_p; the largest TR label is selected.
% Use m-VCF HB/HB and epsilon=0.1 to build prob2 via ode_TR2QuasiP.
% This script currently stops after problem construction: it has no coco call,
% does not compute a new continuation branch, and does not itself draw Fig. 14.
% prob2 remains in the workspace for initialization inspection or further use.
% --------------------------------------------------------------------------
% Read a prerequisite saved run from the current COCO data directory.
bd = coco_bd_read('qp_d1');
% Extract saved Neimark-Sacker/torus event labels from the periodic branch.
lab  = coco_bd_labs(bd, 'TR');
lab  = max(lab);
prob2 = coco_prob();
% Reconstruct displacement/velocity output for the selected coordinates.
% Nperiod/Nsampling/Omt determine the requested time-domain sampling.
amp_opts = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2,3],'Nperiod',1,'Nsampling',2^8);
% For d=1 evaluate Floquet stability; for d=2 evaluate Lyapunov stability.
% S/L1/N_Ly control sampling, integration resolution, and iteration length.
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',500,'Omt',[1,2]);
prob2 = coco_set(prob2,'cont','PtMX',3);
prob2 = coco_set(prob2,'cont','NSV',1);
prob2 = coco_set(prob2,'cont','h_min',1e-4,'h_max',0.01);
prob2 = coco_set(prob2,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
% Choose torus parameterization, numerical method, basis orders, and sampling grids.
% MP/mVCF uses HB/FD/CO; MT/SHT uses HT/MF and the chosen ODE integrator.
% P stores period factors; S_wave controls denser output reconstruction.
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'           , ...     
    'type',    { {'HB','HB'} }, ...  
    'k',       { {9,7} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^6,2^6} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
% Predict a 2-torus using the critical multiplier/eigenvector and epsilon=0.1.
% Model callbacks and ODE settings are inherited; the response dimension becomes d=2.
prob2 = ode_TR2QuasiP(prob2, '', 'qp_d1',  lab , 0.1);