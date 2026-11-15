% FIG_5_2_3  Initial quasi-periodic guesses at the two Airfoil TR points.
% QuasiP article, Section 5.2, Fig. 14.
% Run demo_5_2_p first, then run this script from the Airfoil directory.
% Airfoil_TR_predict is the example-local copy of the prediction routines.
% The predicted tori are initial guesses, not corrected solutions.
% This script follows the figure/plot3/plot_setting_paper style of the examples.
%% ====================== TR1 ======================
bd = coco_bd_read('qp_d1');
lab = coco_bd_labs(bd, 'TR');
lab = min(lab);
lab_all = coco_bd_col(bd, 'LAB');
idx_TR = cellfun(@(x) ~isempty(x) && x == lab, lab_all);
bd_TR = bd([1, find(idx_TR)'+1], :);

Psi  = coco_bd_col(bd_TR, 'Psi');
Om1  = coco_bd_col(bd_TR, 'Om');
Q_TR = coco_bd_col(bd_TR, 'Q');
x_po = coco_bd_col(bd_TR, 'HT_x');
u_po = coco_bd_col(bd_TR, 'HT_u');
y_po = [x_po, u_po];
[~, data] = QuasiP_read_solution('', 'qp_d1', lab);
P1 = data.PHI.P{1};
% Same sampling order as the original 256-by-256 predictor grid.
S = min(size(y_po,1), size(Psi,3)-1);
Y = S/P1;

figure(1)
plot3([y_po(:,1);y_po(1,1)], [y_po(:,4);y_po(1,4)], ...
      [y_po(:,2);y_po(1,2)], 'k--', 'LineWidth',1.5)
hold on

% Predict and display the torus for epsilon = 0.05.
[Xk_MF0, OM0, y_qp] = Airfoil_TR_predict(Psi, Om1, P1, 0.05, y_po);
y_qp = reshape(y_qp, Y, S, []);
y_qp = y_qp(1:8:end, 1:8:end, :);
y_qp = reshape(y_qp, [], size(y_po,2));
plot3(y_qp(:,1), y_qp(:,4), y_qp(:,2), 'g.')

% Predict and display the torus for epsilon = 0.10.
[Xk_MF0, OM0, y_qp] = Airfoil_TR_predict(Psi, Om1, P1, 0.10, y_po);
y_qp = reshape(y_qp, Y, S, []);
y_qp = y_qp(1:8:end, 1:8:end, :);
y_qp = reshape(y_qp, [], size(y_po,2));
plot3(y_qp(:,1), y_qp(:,4), y_qp(:,2), 'r.')

% Predict and display the torus for epsilon = 0.20.
[Xk_MF0, OM0, y_qp] = Airfoil_TR_predict(Psi, Om1, P1, 0.20, y_po);
y_qp = reshape(y_qp, Y, S, []);
y_qp = y_qp(1:8:end, 1:8:end, :);
y_qp = reshape(y_qp, [], size(y_po,2));
plot3(y_qp(:,1), y_qp(:,4), y_qp(:,2), 'b.')

xlabel('$q_1$', 'Interpreter','latex');
ylabel('$\dot{q}_1$', 'Interpreter','latex');
zlabel('$q_2$', 'Interpreter','latex');
legend('PO at TR', 'QP-$\varepsilon=0.05$', ...
       'QP-$\varepsilon=0.10$', 'QP-$\varepsilon=0.20$', 'Location','best');
grid on;
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',3)
legend boxon;
hold off

%% ====================== TR2 ======================
bd = coco_bd_read('qp_d1');
lab = coco_bd_labs(bd, 'TR');
lab = max(lab);
lab_all = coco_bd_col(bd, 'LAB');
idx_TR = cellfun(@(x) ~isempty(x) && x == lab, lab_all);
bd_TR = bd([1, find(idx_TR)'+1], :);

Psi  = coco_bd_col(bd_TR, 'Psi');
Om1  = coco_bd_col(bd_TR, 'Om');
Q_TR = coco_bd_col(bd_TR, 'Q');
x_po = coco_bd_col(bd_TR, 'HT_x');
u_po = coco_bd_col(bd_TR, 'HT_u');
y_po = [x_po, u_po];
[~, data] = QuasiP_read_solution('', 'qp_d1', lab);
P1 = data.PHI.P{1};
% Same sampling order as the original 256-by-256 predictor grid.
S = min(size(y_po,1), size(Psi,3)-1);
Y = S/P1;

figure(2)
plot3([y_po(:,1);y_po(1,1)], [y_po(:,4);y_po(1,4)], ...
      [y_po(:,2);y_po(1,2)], 'k--', 'LineWidth',1.5)
hold on

% Predict and display the torus for epsilon = 0.05.
[Xk_MF0, OM0, y_qp] = Airfoil_TR_predict(Psi, Om1, P1, 0.05, y_po);
y_qp = reshape(y_qp, Y, S, []);
y_qp = y_qp(1:8:end, 1:8:end, :);
y_qp = reshape(y_qp, [], size(y_po,2));
plot3(y_qp(:,1), y_qp(:,4), y_qp(:,2), 'g.')

% Predict and display the torus for epsilon = 0.10.
[Xk_MF0, OM0, y_qp] = Airfoil_TR_predict(Psi, Om1, P1, 0.10, y_po);
y_qp = reshape(y_qp, Y, S, []);
y_qp = y_qp(1:8:end, 1:8:end, :);
y_qp = reshape(y_qp, [], size(y_po,2));
plot3(y_qp(:,1), y_qp(:,4), y_qp(:,2), 'r.')

% Predict and display the torus for epsilon = 0.20.
[Xk_MF0, OM0, y_qp] = Airfoil_TR_predict(Psi, Om1, P1, 0.20, y_po);
y_qp = reshape(y_qp, Y, S, []);
y_qp = y_qp(1:8:end, 1:8:end, :);
y_qp = reshape(y_qp, [], size(y_po,2));
plot3(y_qp(:,1), y_qp(:,4), y_qp(:,2), 'b.')

xlabel('$q_1$', 'Interpreter','latex');
ylabel('$\dot{q}_1$', 'Interpreter','latex');
zlabel('$q_2$', 'Interpreter','latex');
legend('PO at TR', 'QP-$\varepsilon=0.05$', ...
       'QP-$\varepsilon=0.10$', 'QP-$\varepsilon=0.20$', 'Location','best');
grid on;
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',3)
legend boxon;
hold off
