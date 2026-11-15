% DEMO_5_1_4_BP_A  Continue embedded Tor11 in the P=(2,2) Tor22 problem.
% Independent experiment: original Duffing demos and runs are unchanged.
% Start from the saved Tor11 Fourier data, with no time integration.
% Retain the Tor22 HB orders (7,7), grid (32,32), and period factors (2,2).
% BP detection uses the algebraic Jacobian; expensive Lyapunov processing
% is disabled during this scan and is not used to label stability.
clearvars;
close all;
clc;
example_dir = fileparts(mfilename('fullpath'));
duffing_dir = fileparts(example_dir);
addpath(duffing_dir);
cd(example_dir);
source_run = fullfile(duffing_dir,'data','qp514_tor11');
bd = coco_bd_read(source_run);
labs = coco_bd_labs(bd);
lab_cells = coco_bd_col(bd,'LAB',false);
% Use the first saved point, at the low-frequency stable side.
source_lab = labs(1);
row = find(cellfun(@(v) isequal(v,source_lab),lab_cells));
bd0 = bd([1,row+1],:);
MF_f = coco_bd_col(bd0,'MF_f');
Xk_MF.k_wave = MF_f(:,3:4); % Physical wave numbers: do not multiply by two.
Xk_MF.Xc = coco_bd_col(bd0,'MF_x');
[sol0,data0] = QuasiP_read_solution('',source_run,source_lab);
OM0 = sol0.u(data0.OM_idx);
p0 = sol0.u(data0.p_idx);

%% ====================== Tor11 in P=(2,2) ======================
prob = coco_prob();
prob = coco_set(prob,'cont','PtMX',500,'NSV',10,'BP',true, ...
    'h0',0.01,'h_min',1e-3,'h_max',0.08,'bi_direct',false);
prob = coco_set(prob,'corr','ItMX',30);
ode_opts = struct('n',1,'d',2,'e',2,'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',false,'isht',false, ...
    'outdof',1,'Nperiod',1,'Nsampling',256);
sta_opts = struct('issta',false);
PHI_INITIAL = struct('param','MP','method','mVCF','isvar',true, ...
    'type',{{'HB','HB'}},'k',{{7,7}},'P',{{2,2}}, ...
    'S',{{32,32}},'S_wave',{{256,256}});
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts, ...
    'sta_opts',sta_opts,'PHI_INITIAL',PHI_INITIAL);
funcs = {@Duffing1_F,@Duffing1_DFDX,@Duffing1_DFDP,@Duffing1_L};
prob = ode_isol2QuasiP(prob,'',funcs{:},Xk_MF,OM0,data0.pnames_np,p0);
[data,uidx] = coco_get_func_data(prob,'QuasiP','data','uidx');
prob = coco_add_func(prob,'KOM_zero',@frequency_ratio,[],'zero', ...
    'uidx',uidx(data.kOM_idx),'F+dF');
coco(prob,'tor11_P22',[],1,{'k_OM1','k_OM2','zeta','gamma','f1','f2'},[1.8,3]);

bd = coco_bd_read('tor11_P22');
bp_labs = coco_bd_labs(bd,'BP');
fprintf('Detected BP labels: %s\n',mat2str(bp_labs));


function [data,y,J] = frequency_ratio(prob,data,u)
y = u(1)/sqrt(2)-u(2);
J = [1/sqrt(2),-1];
end
