% FIG_5_2_1  Plot periodic/quasi-periodic response amplitudes versus Q.
% Reference: Section 5.2, Fig. 12; coordinates are plunge, wing pitch, store pitch.
% Requires qp_d1, qp_d2_mvf_a, qp_d2_mvf_b, and qp_d2_sht_a.
% Run demo_5_2_p, demo_5_2_qp_mVCF, and demo_5_2_qp_SHT first.
% Periodic-branch TR/FP markers and stable/unstable segments are read from bd.
% SHT points are subsampled for display; use full branch data for diagnostics.
% The plot is a comparison of saved results, not a new numerical solve.
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


figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'k:','LineWidth',1.5}; % unstable
ST2{2} = {'k-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_p,max_p(:,1),Stab_p,ST2,'PO-unstable','PO-stable')
h_TR = plot(Q_p(idx_TR),max_p(idx_TR,1),'c^','MarkerSize',8);
h_FP = plot(Q_p(idx_FP),max_p(idx_FP,1),'go','MarkerSize',8);
legend('PO-unstable','PO-stable','PO-TR','PO-FP','Location','best');
xlabel('Q');
ylabel('$\Vert q_1 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'k:','LineWidth',1.5}; % unstable
ST2{2} = {'k-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_p,max_p(:,2),Stab_p,ST2,'PO-unstable','PO-stable')
h_TR = plot(Q_p(idx_TR),max_p(idx_TR,2),'c^','MarkerSize',8);
h_FP = plot(Q_p(idx_FP),max_p(idx_FP,2),'go','MarkerSize',8);
legend('PO-unstable','PO-stable','PO-TR','PO-FP','Location','best');
xlabel('Q');
ylabel('$\Vert q_2 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(3)
hold on
ST2 = cell(2,1);
ST2{1} = {'k:','LineWidth',1.5}; % unstable
ST2{2} = {'k-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_p,max_p(:,3),Stab_p,ST2,'PO-unstable','PO-stable')
h_TR = plot(Q_p(idx_TR),max_p(idx_TR,3),'c^','MarkerSize',8);
h_FP = plot(Q_p(idx_FP),max_p(idx_FP,3),'go','MarkerSize',8);
legend('PO-unstable','PO-stable','PO-TR','PO-FP','Location','best');
xlabel('Q');
ylabel('$\Vert q_3 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


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
plot_stab_lines4(Q_qp,max_qp(:,1),Stab_qp(:,1),ST2,'QP-unstable','QP-stable')
xlabel('Q');
ylabel('$\Vert q_1 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()

figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_qp,max_qp(:,2),Stab_qp(:,1),ST2,'QP-unstable','QP-stable')
xlabel('Q');
ylabel('$\Vert q_2 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()

figure(3)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_qp,max_qp(:,3),Stab_qp(:,1),ST2,'QP-unstable','QP-stable')
xlabel('Q');
ylabel('$\Vert q_3 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
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
plot_stab_lines4(Q_qp,max_qp(:,1),Stab_qp(:,1),ST2)
xlabel('Q');
ylabel('$\Vert q_1 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_qp,max_qp(:,2),Stab_qp(:,1),ST2)
xlabel('Q');
ylabel('$\Vert q_2 \Vert_{\infty}$','Interpreter','latex');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()

figure(3)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Q_qp,max_qp(:,3),Stab_qp(:,1),ST2)
xlabel('Q');
ylabel('$\Vert q_3 \Vert_{\infty}$','Interpreter','latex');
grid on;
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
[Q_qpSample,max_qpSample,Stab_qpSample,~] = ...
    sample_curve_points(Q_qp,max_qp(:,1),Stab_qp,22,0.1);
plot(Q_qpSample, max_qpSample,...
    'bo','MarkerSize',6);
xlabel('Q');
ylabel('$\Vert q_1 \Vert_{\infty}$','Interpreter','latex');
legend('PO-unstable','PO-stable','PO-TR','PO-FP','QP-stable m-VCF','QP-stable FSE-SHT','Location','best');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
[Q_qpSample,max_qpSample,Stab_qpSample,~] = ...
    sample_curve_points(Q_qp,max_qp(:,2),Stab_qp,22,0.1);
plot(Q_qpSample, max_qpSample,...
    'bo','MarkerSize',6);
xlabel('Q');
ylabel('$\Vert q_2 \Vert_{\infty}$','Interpreter','latex');
legend('PO-unstable','PO-stable','PO-TR','PO-FP','QP-stable m-VCF','QP-stable FSE-SHT','Location','best');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()

figure(3)
hold on
[Q_qpSample,max_qpSample,Stab_qpSample,~] = ...
    sample_curve_points(Q_qp,max_qp(:,3),Stab_qp,22,0.1);
plot(Q_qpSample, max_qpSample,...
    'bo','MarkerSize',6);
xlabel('Q');
ylabel('$\Vert q_3 \Vert_{\infty}$','Interpreter','latex');
legend('PO-unstable','PO-stable','PO-TR','PO-FP','QP-stable m-VCF','QP-stable FSE-SHT','Location','best');
grid on;
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()
