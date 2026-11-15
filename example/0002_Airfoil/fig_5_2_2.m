% FIG_5_2_2  Compare periodic and quasi-periodic base frequencies.
% Reference: Section 5.2, Fig. 13 (Q-frequency and frequency-frequency views).
% Requires qp_d1, qp_d2_mvf_a, qp_d2_mvf_b, and qp_d2_sht_a.
% The periodic branch supplies one unknown frequency; a 2-torus supplies two.
% The interval selection assumes at least two saved TR points.
% The selected imaginary Floquet exponent is used as a frequency predictor.
% Its raw array index is not guaranteed to track the same mode along a branch;
% verify eigenvalue continuity before interpreting that curve physically.
% This script plots saved data and does not generate a new continuation run.
% --------------------------------------------------------------------------
% Read a prerequisite saved run from the current COCO data directory.
bd1 = coco_bd_read('qp_d1');
max_p  = coco_bd_col(bd1, 'Max_x');
Om_p   = coco_bd_col(bd1, 'Om');
Q_p    = coco_bd_col(bd1, 'Q');
Ly_p   = coco_bd_col(bd1,'Ly');
Stab_p = coco_bd_col(bd1, 'Stab');
lab    = coco_bd_col(bd1,'LAB');
lab_vec = zeros(length(lab),1);
for ii = 1:length(lab)
    if ~isempty(lab{ii})
        lab_vec(ii) = lab{ii};
    else
        lab_vec(ii) = 0;
    end
end
% Extract saved Neimark-Sacker/torus event labels from the periodic branch.
labTR  = coco_bd_labs(bd1, 'TR');
labFP  = coco_bd_labs(bd1, 'FP');
idx_TR = ismember(lab_vec, labTR);
idx_FP = ismember(lab_vec, labFP);
pos_TR = find(idx_TR);  
num    = pos_TR(1):1:pos_TR(2);


figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'k:','LineWidth',1.5}; % unstable
ST2{2} = {'k-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_p,Om_p,Stab_p,ST2,'PO-unstable','PO-stable')
h_TR = plot(Q_p(idx_TR),Om_p(idx_TR),'c^','MarkerSize',8);
h_FP = plot(Q_p(idx_FP),Om_p(idx_FP),'go','MarkerSize',8);
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
Ly_p = reshape(Ly_p(:,2,num),6,[]);
plot(Om_p(num),imag(Ly_p(1,:)),'k:');
legend('PO-unstable');
hold on


bd1 = coco_bd_read('qp_d2_mvf_a');
max_qp  = coco_bd_col(bd1, 'Max_x');
Om_qp   = coco_bd_col(bd1, 'Om');
Q_qp    = coco_bd_col(bd1, 'Q');
Ly_qp   = coco_bd_col(bd1,'Ly');
Stab_qp = coco_bd_col(bd1, 'Stab');

figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_qp,Om_qp(:,1),Stab_qp(:,1),ST2,'QP-unstable','QP-stable')
plot_stab_lines4(Q_qp,Om_qp(:,2),Stab_qp(:,1),ST2)
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_qp(:,1),Om_qp(:,2),Stab_qp(:,1),ST2,'QP-unstable','QP-stable')% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


bd1 = coco_bd_read('qp_d2_mvf_b');
max_qp  = coco_bd_col(bd1, 'Max_x');
Om_qp   = coco_bd_col(bd1, 'Om');
Q_qp    = coco_bd_col(bd1, 'Q');
Ly_qp   = coco_bd_col(bd1,'Ly');
Stab_qp = coco_bd_col(bd1, 'Stab');

figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_qp,Om_qp(:,1),Stab_qp(:,1),ST2)
plot_stab_lines4(Q_qp,Om_qp(:,2),Stab_qp(:,1),ST2)
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_qp(:,1),Om_qp(:,2),Stab_qp(:,1),ST2)% xlim([1.8 5])
xlabel('$\omega_1$','Interpreter','latex');
ylabel('$\omega_2$','Interpreter','latex');
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()



bd1 = coco_bd_read('qp_d2_sht_a');
max_qp  = coco_bd_col(bd1, 'Max_x');
Om_qp   = coco_bd_col(bd1, 'Om');
Q_qp    = coco_bd_col(bd1, 'Q');
Ly_qp   = coco_bd_col(bd1,'Ly');
Stab_qp = coco_bd_col(bd1, 'Stab');

figure(1)
hold on
[Q_qpSample,Om_qpSample,Stab_qpSample,~] = ...
    sample_curve_points(Q_qp,Om_qp(:,1),Stab_qp,22,0.1);
plot(Q_qpSample, Om_qpSample,...
    'bo','MarkerSize',6);
[Q_qpSample,Om_qpSample,Stab_qpSample,~] = ...
    sample_curve_points(Q_qp,Om_qp(:,2),Stab_qp,22,0.1);
plot(Q_qpSample, Om_qpSample,...
    'bo','MarkerSize',6);
grid on;
xlabel('Q');
ylabel('$\omega_i$','Interpreter','latex');
legend('PO-stable','PO-unstable','PO-TR','PO-FP','QP-stable m-VCF','QP-stable FSE-SHT','Location','best');
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
[Q_qpSample,Om_qp1Sample,Om_qp2Sample,~] = ...
    sample_curve_points(Q_qp,Om_qp(:,1),Om_qp(:,2),22,0.1);
plot(Om_qp1Sample,Om_qp2Sample,...
    'bo','MarkerSize',6);
grid on;
legend('PO-unstable','QP-stable m-VCF','QP-stable FSE-SHT');
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()

