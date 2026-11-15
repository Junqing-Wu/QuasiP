% DEMO_5_2_QP_MVCF  Initialize and continue 2-tori using m-VCF.
% Reference: Section 5.2, Figs. 12-13; requires qp_d1 from demo_5_2_p.
% Read the smallest saved TR label, intended to select the lower TR near Q=7.0762.
% Verify its actual Q: label order is a storage convention, not a physical criterion.
% ode_TR2QuasiP inherits the model and changes the response dimension to d=2.
% Since e=0, both base frequencies are unknown and two OPC rows are needed.
% Use epsilon=0.1 for the critical-eigenvector perturbation of the periodic orbit.
% Two runs start from the selected TR: qp_d2_mvf_a and qp_d2_mvf_b.
% The first has a small point budget near the lower endpoint; the second
% continues over a wider Q interval. Inspect actual coverage in the logs.
% Next: run the SHT comparison and fig_5_2_1 / fig_5_2_2.
% --------------------------------------------------------------------------
%% continuation of QP from TR bifurcation
% Read a prerequisite saved run from the current COCO data directory.
bd = coco_bd_read('qp_d1');
% Extract saved Neimark-Sacker/torus event labels from the periodic branch.
lab  = coco_bd_labs(bd, 'TR');
lab  = min(lab);
prob2 = coco_prob();
% Reconstruct displacement/velocity output for the selected coordinates.
% Nperiod/Nsampling/Omt determine the requested time-domain sampling.
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2,3],'Nperiod',1,'Nsampling',2^8);
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
% Continue the one-dimensional branch in Q; inspect actual coverage and stopping reasons.
coco(prob2, 'qp_d2_mvf_a', [], 1, { 'Q' 'u_om1' 'u_om2'}, [7.07 7.07615045733112]);

bd = coco_bd_read('qp_d1');
lab  = coco_bd_labs(bd, 'TR');
lab  = min(lab);
prob2 = coco_prob();
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2,3],'Nperiod',1,'Nsampling',2^8);
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',500,'Omt',[1,2]);
prob2 = coco_set(prob2,'cont','PtMX',59);
prob2 = coco_set(prob2,'cont','NSV',1);
prob2 = coco_set(prob2,'cont','h_min',1e-3,'h_max',0.1);
prob2 = coco_set(prob2,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'           , ...     
    'type',    { {'HB','HB'} }, ...  
    'k',       { {9,7} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^6,2^6} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
prob2 = ode_TR2QuasiP(prob2, '', 'qp_d1',  lab , 0.1);
% Continue the one-dimensional branch in Q; inspect actual coverage and stopping reasons.
coco(prob2, 'qp_d2_mvf_b', [], 1, { 'Q' 'u_om1' 'u_om2'}, [7.07615045733112 10]);