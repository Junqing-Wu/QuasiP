% FIGURE_RS_FRC  Periodic/quasi-periodic frequency-response curves with stability.
% Run demo_TR2QP first, then execute this script from the RSrubbing directory.
% qp1 contains periodic orbits; qp2 contains quasi-periodic tori.
% Amplitudes are max(abs(x)) and max(abs(y)) over the reconstructed orbit/torus.
% Black: PO; red: QP. Solid: stable; dashed: unstable; dotted: undetermined.
% Stability is read from Stab; this script does not recompute exponents.
% Preserve continuation order: sorting by frequency would corrupt folded branches.

%% ====================== Read the periodic branch ======================
bd_p = coco_bd_read('qp1');
Om_p = coco_bd_col(bd_p,'k_OM1');
max_p = coco_bd_col(bd_p,'Max_x');
Stab_p = coco_bd_col(bd_p,'Stab');
Om_p = Om_p(:);
S_p = frc_stability_state(Stab_p,numel(Om_p));

%% ====================== Read the quasi-periodic branch ======================
bd_qp = coco_bd_read('qp2');
Om_qp = coco_bd_col(bd_qp,'k_OM1');
max_qp = coco_bd_col(bd_qp,'Max_x');
Stab_qp = coco_bd_col(bd_qp,'Stab');
Om_qp = Om_qp(:);
S_qp = frc_stability_state(Stab_qp,numel(Om_qp));
assert(size(max_p,1)==numel(Om_p) && size(max_qp,1)==numel(Om_qp), ...
    'Amplitude arrays do not match their continuation branches.');
assert(size(max_p,2)>=2 && size(max_qp,2)>=2, ...
    'Both displacement amplitudes are required. Rerun with outdof=[1,2].');

%% ====================== x-displacement FRC ======================
figure(1)
hold on
ST2 = cell(3,1);
ST2{1} = {'k--','LineWidth',1.5}; % unstable PO
ST2{2} = {'k-','LineWidth',1.5};  % stable PO
ST2{3} = {'k:','LineWidth',1.5};  % undetermined PO
plot_stab_lines4(Om_p,max_p(:,1),S_p,ST2,'PO-unstable','PO-stable','PO-undetermined')
ST2{1} = {'r--','LineWidth',1.5}; % unstable QP
ST2{2} = {'r-','LineWidth',1.5};  % stable QP
ST2{3} = {'r:','LineWidth',1.5};  % undetermined QP
plot_stab_lines4(Om_qp,max_qp(:,1),S_qp,ST2,'QP-unstable','QP-stable','QP-undetermined')
xlabel('$\omega_1$','Interpreter','latex');
ylabel('$\Vert x\Vert_{\infty}$','Interpreter','latex');
grid on;
hold off
plot_setting_paper()
legend boxon;

%% ====================== y-displacement FRC ======================
figure(2)
hold on
ST2 = cell(3,1);
ST2{1} = {'k--','LineWidth',1.5};
ST2{2} = {'k-','LineWidth',1.5};
ST2{3} = {'k:','LineWidth',1.5};
plot_stab_lines4(Om_p,max_p(:,2),S_p,ST2,'PO-unstable','PO-stable','PO-undetermined')
ST2{1} = {'r--','LineWidth',1.5};
ST2{2} = {'r-','LineWidth',1.5};
ST2{3} = {'r:','LineWidth',1.5};
plot_stab_lines4(Om_qp,max_qp(:,2),S_qp,ST2,'QP-unstable','QP-stable','QP-undetermined')
xlabel('$\omega_1$','Interpreter','latex');
ylabel('$\Vert y\Vert_{\infty}$','Interpreter','latex');
grid on;
hold off
plot_setting_paper()
legend boxon;

function S = frc_stability_state(Stab,npoints)
% plot_stab_lines4 uses S+1 as the style index: 0=unstable, 1=stable, 2=unknown.
% If several Omt estimates are stored, require agreement for a binary label.
if isvector(Stab) && numel(Stab)==npoints
    Stab = reshape(Stab,npoints,1);
end
assert(size(Stab,1)==npoints,'Stability data do not match the branch point count.');
S = 2*ones(npoints,1);
S(all(isfinite(Stab) & Stab==0,2)) = 0;
S(all(isfinite(Stab) & Stab==1,2)) = 1;
end
