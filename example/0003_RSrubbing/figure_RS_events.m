% FIGURE_RS_EVENTS  XY orbits at prescribed rotor speeds.
% Run demo_TR2QP first to save qp1 (periodic) and qp2 (quasi-periodic),
% including the Om1=0.30/0.35/0.40/0.45 event solutions.
% Run this script from the RSrubbing example directory after COCO startup
% and QuasiP install. All event labels are plotted, including multiple
% intersections of a folded branch with the same prescribed frequency.
% No direct integration is overlaid: an unstable continued solution need
% not be the attractor reached by direct numerical integration.

event_names = {'Om1=0.30','Om1=0.35','Om1=0.40','Om1=0.45'};
event_values = [0.30,0.35,0.40,0.45];
run_names = {'qp1','qp2'};
solution_names = {'PO','QP'};
plot_stride = 4;  % Thin the displayed torus samples, not the computed solution.
% Track the requested classes separately; never filter out unstable solutions.
class_counts = zeros(numel(event_values),3); % unstable PO, stable QP, unstable QP

% One XY figure per prescribed frequency. Both runs and all event labels
% contribute to these shared figures. No spectra or velocity portraits are drawn.
for ievent = 1:numel(event_values)
    figure(100+ievent)
    clf
    hold on
end

for irun = 1:numel(run_names)
    run_name = run_names{irun};
    % Read saved branch labels; a run must exist before plotting its events.
    bd = coco_bd_read(run_name);
    for ievent = 1:numel(event_names)
        labs = coco_bd_labs(bd,event_names{ievent});
        if isempty(labs)
            warning('figure_RS_events:MissingEvent', ...
                '%s has no %s solution. Check actual coverage and rerun with events enabled.', ...
                run_name,event_names{ievent});
            continue
        end
        for ilab = 1:numel(labs)
            lab = labs(ilab);
            all_labs = coco_bd_col(bd,'LAB',false);
            row = find(cellfun(@(x) ~isempty(x) && isequal(x,lab),all_labs));
            assert(isscalar(row),'Expected one branch-data row per LAB.');
            event_bd = bd([1,row+1],:);
            stab = coco_bd_col(event_bd,'Stab');
            if ~isempty(stab) && all(isfinite(stab(:))) && all(stab(:)==1)
                stability_name = 'stable';
            elseif ~isempty(stab) && all(isfinite(stab(:))) && all(stab(:)==0)
                stability_name = 'unstable';
            else
                stability_name = 'undetermined';
            end
            [sol,data] = QuasiP_read_solution('',run_name,lab);
            C = sol.u(data.C_idx);
            Om = sol.u(data.OM_idx);
            p = sol.u(data.p_idx);
            assert(abs(Om(1)-event_values(ievent))<1e-5, ...
                'Saved event frequency differs from its requested value.');
            assert(data.ode_opts.n==2 && data.ODE==2 && strcmp(data.PHI.param,'MP'), ...
                'This plot expects the two-coordinate second-order MP rotor model.');
            % Reconstruct both displacement coordinates from the saved
            % coefficients. This also supports older runs whose outdof was 1.
            U_prod = prod(cell2mat(data.PHI.U_omod));
            C_state = reshape(C,U_prod,data.ode_opts.n)';
            HT_x = zeros(prod(cell2mat(data.PHI.pwS_omod)),2);
            for idof = 1:2
                HT_x(:,idof) = MP_mVCF_U2S_x(C_state(idof,:),data.PHI,U_prod);
            end
            curve_name = sprintf('%s-%s',solution_names{irun},stability_name);
            tag = sprintf('%s: k_{OM1}=%.2f, LAB=%d',curve_name,Om(1),lab);
            if data.ode_opts.d==1
                if strcmp(stability_name,'unstable')
                    phase_style='k--'; class_counts(ievent,1)=class_counts(ievent,1)+1;
                else
                    phase_style='k-';
                end
            elseif strcmp(stability_name,'stable')
                phase_style='r.'; class_counts(ievent,2)=class_counts(ievent,2)+1;
            elseif strcmp(stability_name,'unstable')
                phase_style='b.'; class_counts(ievent,3)=class_counts(ievent,3)+1;
            else
                phase_style='m.';
            end

            %% ====================== Rotor-plane orbit ======================
            figure(100+ievent)
            hold on
            orbit_label = sprintf('%s, LAB=%d',curve_name,lab);
            if data.ode_opts.d==1
                plot([HT_x(:,1);HT_x(1,1)],[HT_x(:,2);HT_x(1,2)], ...
                    phase_style,'LineWidth',1.5,'DisplayName',orbit_label)
            else
                plot(HT_x(1:plot_stride:end,1),HT_x(1:plot_stride:end,2), ...
                    phase_style,'MarkerSize',5,'DisplayName',orbit_label)
            end

            fprintf('%s: plotted XY orbit (mu=%.6g, OM=%.6g).\n',tag,p(1),p(2));
        end
    end
end

%% ====================== Format the shared XY figures ======================
for ievent = 1:numel(event_values)
    figure(100+ievent)
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    title(sprintf('k_{OM1}=%.2f',event_values(ievent)));
    axis equal;
    grid on;
    if ~isempty(findobj(gca,'Type','line'))
        legend('show','Location','best');
    end
    hold off
    plot_setting_paper()
    set(findall(gcf,'type','line'),'linewidth',1)
    legend boxon;
end

% Missing classes indicate incomplete branch coverage or absent intersections,
% not permission to relabel a solution. Check the actual continuation runs.
class_names = {'unstable PO','stable QP','unstable QP'};
for ievent=1:numel(event_values)
    for iclass=1:numel(class_names)
        if class_counts(ievent,iclass)==0
            warning('figure_RS_events:MissingClass', ...
                'k_OM1=%.2f: no saved %s event solution was found.', ...
                event_values(ievent),class_names{iclass});
        end
    end
end

