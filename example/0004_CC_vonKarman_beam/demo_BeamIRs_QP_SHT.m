% DEMO_BEAMIRS_QP_SHT  Continue the FSE shooting torus from the lower TR.
% Reference: "Parallelized computation of quasi-periodic solutions for finite
% element problems: A Fourier series expansion-based shooting method".
% Section 4.2, Fig. 16. Run demo_BeamIRs_QP_SHT_seed first.
clearvars;
close all;
clc;

%% ====================== Continue the corrected SHT torus ======================
seedRun = 'beamIRs_qp_sht_seed';
seedLabels = coco_bd_labs(seedRun);
assert(numel(seedLabels)==1,'Run demo_BeamIRs_QP_SHT_seed first.');
targets = [33.859,33.931,34.034];
prob = coco_prob();
prob = coco_set(prob,'cont','PtMX',[0,100],'NSV',1,'BP',false,'FP',false, ...
    'h0',0.55,'h_max',2.0,'h_min',1e-5,'al_max',18,'TOL',1e-8);
prob = coco_set(prob,'corr','ItMX',16,'ItNW',1,'TOL',2e-7);
[sol,data] = QuasiP_read_solution('',seedRun,seedLabels(1));
% Keep the full state history for amplitudes and only terminal sensitivities.
data.PHI.endpoint_sensitivities = true;
prob = ode_construct_QuasiP(prob,'QuasiP',data,sol);
prob = coco_add_event(prob,'SHT1','k_OM1',targets(1));
prob = coco_add_event(prob,'SHT2','k_OM1',targets(2));
prob = coco_add_event(prob,'SHT3','k_OM1',targets(3));
bd = coco(prob,'beamIRs_qp_sht',[],1,{'k_OM1','u_om2','force_scale'}, ...
    [33.85,34.036]);
