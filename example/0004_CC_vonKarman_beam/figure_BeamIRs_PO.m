% FIGURE_BEAMIRS_PO  Periodic response and Floquet stability.
% Run demo_BeamIRs_PO first. TR is the Neimark-Sacker event in QuasiP.
runName = 'beamIRs_po_stab';
bd = coco_bd_read(runName);
omega = coco_bd_col(bd,'k_OM1');
amp = coco_bd_col(bd,'Max_x');
stab = coco_bd_col(bd,'Stab');
ly = coco_bd_col(bd,'Ly',false);
rho = cellfun(@(a) max(abs(a(:,1))),ly);
labs = coco_bd_col(bd,'LAB',false);
savedLabels = coco_bd_labs(bd);
[~,data] = QuasiP_read_solution('',runName,savedLabels(1));
% Stability processing saves Max_x for all DOFs in model order.
amp = amp(:,data.amp_opts.outdof);
folds = coco_bd_labs(bd,'FP');
tori = coco_bd_labs(bd,'TR');
doubles = coco_bd_labs(bd,'PD');
isfold = cellfun(@(lab) any(ismember(lab,folds)),labs);
istorus = cellfun(@(lab) any(ismember(lab,tori)),labs);
isdouble = cellfun(@(lab) any(ismember(lab,doubles)),labs);

% Preserve continuation order through both folds.
% Use the classification after an event to start the outgoing edge.
edge = stab(1:end-1);
atEvent = isfold | istorus | isdouble;
next = stab(2:end);
edge(atEvent(1:end-1)) = next(atEvent(1:end-1));
edgeBreak = [1,find(diff(edge)~=0)+1,numel(omega)];
styles = {'b--','b-'};
legendNames = {'Unstable PO','Stable PO'};
yLabels = {'$\|w(l/4)\|_\infty$ [mm]','$\|w(l/2)\|_\infty$ [mm]'};

for idof = 1:2
    figure(idof)
    clf
    hold on
    shown = false(1,2);
    for iseg = 1:numel(edgeBreak)-1
        idx = edgeBreak(iseg):edgeBreak(iseg+1);
        style = edge(idx(1))+1;
        h = plot(omega(idx),amp(idx,idof),styles{style},'DisplayName',legendNames{style});
        if shown(style)
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        shown(style) = true;
    end
    if any(isfold)
        plot(omega(isfold),amp(isfold,idof),'ro','MarkerSize',7,'DisplayName','Fold');
    end
    if any(istorus)
        plot(omega(istorus),amp(istorus,idof),'ms','MarkerSize',7,'DisplayName','NS (TR)');
    end
    if any(isdouble)
        plot(omega(isdouble),amp(isdouble,idof),'kd','MarkerSize',7,'DisplayName','PD');
    end
    xlabel('$\omega_1$ [rad/s]','Interpreter','latex');
    ylabel(yLabels{idof},'Interpreter','latex');
    grid on
    legend('show','Location','best');
    plot_setting_paper()
    hold off
end

figure(3)
clf
plot(omega,rho,'k-','DisplayName','Largest modulus')
hold on
yline(1,'k:','DisplayName','Unit circle');
if any(istorus)
    plot(omega(istorus),rho(istorus),'ms','MarkerSize',7,'DisplayName','NS (TR)');
end
xlabel('$\omega_1$ [rad/s]','Interpreter','latex');
ylabel('$\max_j|\mu_j|$','Interpreter','latex');
grid on
legend('show','Location','best');
plot_setting_paper()
hold off
