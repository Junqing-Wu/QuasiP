% DEMO_5_2_QP_SHT  Continue 2-tori using Fourier-series-expansion shooting.
% Reference: Section 5.2, Figs. 12-13; requires qp_d1 from demo_5_2_p.
% Read the largest saved TR label, intended to select the upper TR near Q=9.5035.
% Check the saved Q before interpreting the selected label as the upper endpoint.
% The active configuration below uses MT/SHT in the multi-frequency (MF) domain.
% The model is autonomous (e=0), with two unknown frequencies (d=2).
% The current active integrator is ode15s; the commented block is an inactive
% alternative m-VCF configuration and is not executed.
% Output: qp_d2_sht_a. Compare its amplitudes and frequencies with m-VCF
% using fig_5_2_1 and fig_5_2_2 after all prerequisite runs exist.
% Near TR, verify that the solution retains nonzero second-frequency content
% rather than converging back to the periodic orbit, as discussed in Section 5.2.
% --------------------------------------------------------------------------
%% continuation of QP from TR bifurcation
% Read a prerequisite saved run from the current COCO data directory.
bd = coco_bd_read('qp_d1');
% Extract saved Neimark-Sacker/torus event labels from the periodic branch.
lab  = coco_bd_labs(bd, 'TR');
lab  = max(lab);
prob2 = coco_prob();
% Reconstruct displacement/velocity output for the selected coordinates.
% Nperiod/Nsampling/Omt determine the requested time-domain sampling.
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2,3],'Nperiod',1,'Nsampling',2^8);
% For d=1 evaluate Floquet stability; for d=2 evaluate Lyapunov stability.
% S/L1/N_Ly control sampling, integration resolution, and iteration length.
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',500,'Omt',1);
prob2 = coco_set(prob2,'cont','PtMX',42);
prob2 = coco_set(prob2,'cont','NSV',1);
prob2 = coco_set(prob2,'corr','ItMX',50);
prob2 = coco_set(prob2,'cont','h_min',1e-4,'h_max',0.2);
prob2 = coco_set(prob2,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
% Choose torus parameterization, numerical method, basis orders, and sampling grids.
% MP/mVCF uses HB/FD/CO; MT/SHT uses HT/MF and the chosen ODE integrator.
% P stores period factors; S_wave controls denser output reconstruction.
PHI_INITIAL = struct( 'param', 'MT', ...
    'method',     'SHT'           , ...     
    'type',       'MF' , ... % HT or MF
    'ODEsolver',  @ode15s, ...
    'Omt',        1, ...
    'M',          1, ...
    'L1',         2^8,...
    'k',          { {9,5} }, ...        
    'P',          { {1,1} }, ...    
    'S',          { {2^6,2^6} }, ...
    'S_wave',     { {2^8,2^8} });
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
% Predict a 2-torus using the critical multiplier/eigenvector and epsilon=0.1.
% Model callbacks and ODE settings are inherited; the response dimension becomes d=2.
prob2 = ode_TR2QuasiP(prob2, '', 'qp_d1',  lab , 0.1);
% Continue the one-dimensional branch in Q; inspect actual coverage and stopping reasons.
coco(prob2, 'qp_d2_sht_a', [], 1, { 'Q' 'u_om1' 'u_om2'}, [7.076 9.50353970388679]);

% bd = coco_bd_read('qp_d1');
% lab  = coco_bd_labs(bd, 'TR');
% lab  = min(lab);
% prob2 = coco_prob();
% amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2,3],'Nperiod',1,'Nsampling',2^8);
% sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',500,'Omt',[1,2]);
% prob2 = coco_set(prob2,'cont','PtMX',59);
% prob2 = coco_set(prob2,'cont','NSV',1);
% prob2 = coco_set(prob2,'cont','h_min',1e-3,'h_max',0.1);
% prob2 = coco_set(prob2,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
% PHI_INITIAL = struct( 'param', 'MP', ...
%     'method',  'mVCF'           , ...     
%     'type',    { {'HB','FD'} }, ...  
%     'k',       { {9,-2:1:1} }      , ...        
%     'P',       { {1,1} }      , ...    
%     'S',       { {2^5,2^5} }  , ...
%     'S_wave',  { {2^8,2^8} }  , ...
%     'OmOd' ,   [2,1]);
% prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
% prob2 = ode_TR2QuasiP(prob2, '', 'qp_d1',  lab , 0.1);
% coco(prob2, 'qp_d2_mvf_bb', [], 1, { 'Q' 'u_om1' 'u_om2'}, [7.07615045733112 10]);