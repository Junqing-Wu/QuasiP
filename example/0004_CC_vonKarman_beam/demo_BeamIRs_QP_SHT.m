% DEMO_BEAMIRS_QP_SHT  Correct and continue the torus from the lower TR.
% Reference: "Parallelized computation of quasi-periodic solutions for finite
% element problems: A Fourier series expansion-based shooting method".
% Section 4.2, Fig. 16. Requires beamIRs_po_stab from demo_BeamIRs_PO.
clearvars;
close all;
clc;

%% ====================== TR initial guess ======================
poRun = 'beamIRs_po_stab';
bd = coco_bd_read(poRun);
labs = coco_bd_labs(bd,'TR');
assert(numel(labs)>=2,'Compute the periodic branch with TR detection first.');
lab = labs(1);
saved = coco_bd_col(bd,'LAB',false);
row = find(cellfun(@(v) isequal(v,lab),saved));
psi = coco_bd_col(bd([1,row+1],:),'Psi');
[~,data] = QuasiP_read_solution('',poRun,lab);
[V,mu] = eig(psi(:,:,end),'vector');
pair = find(imag(mu)>1e-6);
[~,idx] = min(abs(abs(mu(pair))-1));
v = V(:,pair(idx));
outdof = data.amp_opts.outdof(1);
% Scale the unit-norm Floquet eigenvector to a 0.17 mm trial modulation.
scale = max(abs(squeeze(pagemtimes(psi(outdof,:,:),v))));
epsilon = 0.17/scale;
%% ====================== SHT seed correction ======================
prob = coco_prob();
amp_opts = struct('isamp',true,'isuamp',false,'isht',false, ...
    'outdof',1:data.ode_opts.n,'Nperiod',0,'Nsampling',256);
sta_opts = struct('issta',false);
PHI_INITIAL = struct('param','MT','method','SHT','type','MF', ...
    'ODEsolver',@Newmark,'ODE2slover_opt',struct('ItMX',12,'RelTol',1e-4,'alpha',0), ...
    'Omt',1,'M',1,'L1',1024,'residual_only',true,'endpoint_sensitivities',true,'k',{{7,7}},'P',{{1,1}}, ...
    'S',{{1024,32}},'S_wave',{{128,32}});
par_opts = struct('parallel',true,'ncores',8);
prob = coco_set(prob,'QuasiP','amp_opts',amp_opts, ...
    'sta_opts',sta_opts,'par_opts',par_opts,'PHI_INITIAL',PHI_INITIAL);
prob = coco_set(prob,'cont','PtMX',[0,100],'NSV',1,'BP',false,'FP',false, ...
    'h0',0.02,'h_max',0.12,'h_min',1e-5,'TOL',1e-8);
prob = coco_set(prob,'corr','ItMX',16,'TOL',2e-7);
% Fix a sideband projection for the first correction so the torus cannot
% collapse to the periodic family. The excitation frequency can adjust.
prob = ode_TR2QuasiP(prob,'',poRun,lab,epsilon);
[data,uidx,u0] = coco_get_func_data(prob,'QuasiP','data','uidx','u0');
idx = reshape(data.C_idx,data.PHI.U,[]);
idx = idx(2:3,1:data.ode_opts.n);
coeff = u0(idx(:));
seedFix = struct('direction',coeff'/norm(coeff),'target',norm(coeff));
prob = coco_add_func(prob,'seed_sideband',@fix_sideband,seedFix,'zero', ...
    'uidx',uidx(idx(:)),'F+dF');
seedRun = 'beamIRs_qp_sht_seed';
seed = coco(prob,seedRun,[],0,{'k_OM1','u_om2','force_scale'});
seedLabels = coco_bd_labs(seed);
assert(numel(seedLabels)==1,'The SHT seed correction failed.');
seedWaves = coco_bd_col(seed,'MF_f');
seedCoeffs = coco_bd_col(seed,'MF_x');
assert(norm(seedCoeffs(seedWaves(:,4)~=0,outdof))>1e-3, ...
    'The SHT corrector returned the embedded periodic orbit.');

%% ====================== Continue the corrected SHT torus ======================
targets = [33.859,33.931,34.034];
prob = coco_prob();
prob = coco_set(prob,'cont','PtMX',[0,100],'NSV',1,'BP',false,'FP',false, ...
    'h0',0.55,'h_max',2.0,'h_min',1e-5,'al_max',18,'TOL',1e-8);
prob = coco_set(prob,'corr','ItMX',16,'ItNW',1,'TOL',2e-7);
[sol,data] = QuasiP_read_solution('',seedRun,seedLabels(1));
data.par_opts = par_opts;
% Keep the full state history for amplitudes and only terminal sensitivities.
data.PHI.endpoint_sensitivities = true;
prob = ode_construct_QuasiP(prob,'QuasiP',data,sol);
prob = coco_add_event(prob,'SHT1','k_OM1',targets(1));
prob = coco_add_event(prob,'SHT2','k_OM1',targets(2));
prob = coco_add_event(prob,'SHT3','k_OM1',targets(3));
bd = coco(prob,'beamIRs_qp_sht',[],1,{'k_OM1','u_om2','force_scale'}, ...
    [33.85,34.036]);

function [data,y,J] = fix_sideband(~,data,u)
J = data.direction;
y = J*u-data.target;
end
