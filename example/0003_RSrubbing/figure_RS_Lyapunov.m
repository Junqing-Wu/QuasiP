% FIGURE_RS_LYAPUNOV  Evolution of state stability exponents on the PO/QP branches.
% Run demo_TR2QP first. Execute this script from the RSrubbing directory.
% qp1: Ly stores [Floquet multipliers, Floquet exponents]; plot real(column 2).
% qp2: Ly stores Lyapunov exponents, grouped by sta_opts.Omt.
% The forced rotor has four physical state variables [x,y,xdot,ydot].
% Do not add neutral exponents for excitation phases to these saved state exponents.
% Curves are ranked at each point, not tracked as individual eigenmodes.
% Keep the continuation order to retain folds and separate branch passages.

use_continuation_number = false; % true: point index; false: k_OM1 as in the FRC.

%% ====================== Periodic-orbit exponents ======================
bd_p = coco_bd_read('qp1');
Om_p = coco_bd_col(bd_p,'k_OM1');
Ly_p = coco_bd_col(bd_p,'Ly',false);
Om_p = Om_p(:);
sigma_p = zeros(numel(Ly_p),4);
for ii = 1:numel(Ly_p)
    assert(isequal(size(Ly_p{ii}),[4,2]), ...
        'Expected four Floquet multipliers and four exponents in periodic Ly data.');
    sigma_p(ii,:) = sort(real(Ly_p{ii}(:,2)),'descend')';
end

%% ====================== Quasi-periodic exponents ======================
bd_qp = coco_bd_read('qp2');
Om_qp = coco_bd_col(bd_qp,'k_OM1');
Ly_qp = coco_bd_col(bd_qp,'Ly',false);
Om_qp = Om_qp(:);
labs = coco_bd_labs(bd_qp);
assert(~isempty(labs),'No saved quasi-periodic solution is available.');
[~,data_qp] = QuasiP_read_solution('','qp2',labs(1));
Omt = data_qp.sta_opts.Omt;
sigma_qp = zeros(numel(Ly_qp),4,numel(Omt));
for ii = 1:numel(Ly_qp)
    assert(numel(Ly_qp{ii})==4*numel(Omt), ...
        'The number of saved exponents does not match the Omt directions.');
    values = reshape(real(Ly_qp{ii}),4,[]);
    sigma_qp(ii,:,:) = reshape(sort(values,1,'descend'),1,4,[]);
end

if use_continuation_number
    xp = (1:numel(Om_p))';
    xqp = (1:numel(Om_qp))';
    x_label = 'Continuation point index';
else
    xp = Om_p;
    xqp = Om_qp;
    x_label = '$\omega_1$';
end
colors = lines(4);
line_styles = {'-','--',':','-.'};
% Literal LaTeX labels avoid sprintf interpreting backslash escape sequences.
po_labels = {'$\mathrm{Re}(\lambda_1)$','$\mathrm{Re}(\lambda_2)$', ...
             '$\mathrm{Re}(\lambda_3)$','$\mathrm{Re}(\lambda_4)$'};
qp_labels = {'$\sigma_1$','$\sigma_2$','$\sigma_3$','$\sigma_4$'};

%% ====================== PO spectrum ======================
figure(3)
clf
hold on
for kk = 1:4
    plot(xp,sigma_p(:,kk),'-','Color',colors(kk,:),'LineWidth',1.5, ...
        'DisplayName',po_labels{kk});
end
yline(0,'k:','HandleVisibility','off');
xlabel(x_label,'Interpreter','latex');
ylabel('$\mathrm{Re}(\lambda)$','Interpreter','latex');
title('Periodic solutions');
legend('show','Location','best');
grid on;
hold off
plot_setting_paper()
legend boxon;

%% ====================== QP Lyapunov spectrum ======================
figure(4)
clf
hold on
for idir = 1:numel(Omt)
    style = line_styles{mod(idir-1,numel(line_styles))+1};
    for kk = 1:4
        plot(xqp,sigma_qp(:,kk,idir),style,'Color',colors(kk,:),'LineWidth',1.5, ...
            'DisplayName',[qp_labels{kk},', Omt=',num2str(Omt(idir))]);
    end
end
yline(0,'k:','HandleVisibility','off');
yline(0.001,'k--','DisplayName','QP code threshold');
xlabel(x_label,'Interpreter','latex');
ylabel('$\sigma$','Interpreter','latex');
title('Quasi-periodic solutions');
legend('show','Location','best');
grid on;
hold off
plot_setting_paper()
legend boxon;

%% ====================== Largest exponent comparison ======================
figure(5)
clf
plot(xp,sigma_p(:,1),'k-','LineWidth',1.5,'DisplayName','PO');
hold on
for idir = 1:numel(Omt)
    style = line_styles{mod(idir-1,numel(line_styles))+1};
    plot(xqp,sigma_qp(:,1,idir),style,'Color','r','LineWidth',1.5, ...
        'DisplayName',sprintf('QP, Omt=%d',Omt(idir)));
end
yline(0,'k:','DisplayName','Zero');
yline(0.001,'b--','DisplayName','QP code threshold');
xlabel(x_label,'Interpreter','latex');
ylabel('Largest state exponent');
legend('show','Location','best');
grid on;
hold off
plot_setting_paper()
legend boxon;
