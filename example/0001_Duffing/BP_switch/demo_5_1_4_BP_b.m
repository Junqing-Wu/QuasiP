% DEMO_5_1_4_BP_B  Follow the saved secondary BP vector directly.
% Compare with demo_5_1_4_BP_control, which uses the unmodified default mode.
% The toolbox restart follows the saved secondary BP tangent.
clearvars;
example_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(example_dir));
cd(example_dir);
bd = coco_bd_read('tor11_P22');
bp_labs = coco_bd_labs(bd,'BP');
assert(~isempty(bp_labs),'Run demo_5_1_4_BP_a first.');
%% ====================== Continue from the saved BP points ======================
for ibp = 1:numel(bp_labs)
    prob = coco_prob();
    prob = coco_set(prob,'cont','PtMX',180,'NSV',1, ...
        'h0',0.1,'h_min',1e-4,'h_max',2);
    prob = coco_set(prob,'corr','ItMX',40);
    prob = ode_BP2QuasiP(prob,'','tor11_P22',bp_labs(ibp));
    [data,uidx] = coco_get_func_data(prob,'QuasiP','data','uidx');
    prob = coco_add_func(prob,'KOM_zero',@frequency_ratio,[],'zero', ...
        'uidx',uidx(data.kOM_idx),'F+dF');
    run_name = sprintf('direct_BP%d',bp_labs(ibp));
    coco(prob,run_name,[],1,{'k_OM1','k_OM2','zeta','gamma','f1','f2'},[1.926896,2.830878]);
end

function [data,y,J] = frequency_ratio(prob,data,u)
y = u(1)/sqrt(2)-u(2);
J = [1/sqrt(2),-1];
end
