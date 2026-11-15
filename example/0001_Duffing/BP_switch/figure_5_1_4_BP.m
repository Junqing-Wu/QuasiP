% FIGURE_5_1_4_BP  Check whether BP switching leaves the Tor11 subspace.
% Integer physical wave numbers describe Tor11 even in the P=(2,2) box.
% Noninteger wave numbers are the additional doubled-torus content.
% This test checks branch identity, not Lyapunov stability or resolution convergence.
clearvars;
example_dir = fileparts(mfilename('fullpath'));
cd(example_dir);
bd0 = coco_bd_read('tor11_P22');
bp_labs = coco_bd_labs(bd0,'BP');
all_labs = coco_bd_col(bd0,'LAB',false);
all_om = coco_bd_col(bd0,'k_OM1');
for j=1:numel(bp_labs)
    row = find(cellfun(@(v) isequal(v,bp_labs(j)),all_labs));
    [bp_sol,~] = QuasiP_read_solution('','tor11_P22',bp_labs(j),'BP');
    assert(~isempty(bp_sol.t0),'Missing BP restart vector.');
    fprintf('BP LAB=%d: OM1=%.15g, restart-vector norm=%.6g\n', ...
        bp_labs(j),all_om(row),norm(bp_sol.t0));
end
run_names = {'tor11_P22'};
for j=1:numel(bp_labs)
    for prefix = {'switch','direct'}
        candidate = sprintf('%s_BP%d',prefix{1},bp_labs(j));
        if exist(fullfile('data',candidate,'bd.mat'),'file')
            run_names{end+1} = candidate;
        end
    end
end
figure(1)
clf
hold on
figure(2)
clf
hold on
summary = struct([]);
for irun=1:numel(run_names)
    bd = coco_bd_read(run_names{irun});
    Om = coco_bd_col(bd,'k_OM1');
    amplitude = coco_bd_col(bd,'Max_x');
    MF_f = coco_bd_col(bd,'MF_f',false);
    MF_x = coco_bd_col(bd,'MF_x',false);
    half_norm = zeros(numel(Om),1);
    direction_norm = zeros(numel(Om),2);
    for ii=1:numel(Om)
        k = MF_f{ii}(:,3:4);
        ishalf = abs(k-round(k))>1e-8;
        half_norm(ii) = norm(MF_x{ii}(any(ishalf,2),:),'fro');
        for jj=1:2
            direction_norm(ii,jj) = norm(MF_x{ii}(ishalf(:,jj),:),'fro');
        end
    end
    summary(irun).run = run_names{irun};
    summary(irun).Om = Om;
    summary(irun).half_norm = half_norm;
    summary(irun).direction_norm = direction_norm;
    fprintf('%s: %d points, OM1=[%.12g,%.12g], max half-norm=%.6g, directional=%s\n', ...
        run_names{irun},numel(Om),min(Om),max(Om),max(half_norm), ...
        mat2str(max(direction_norm,[],1),6));
    figure(1)
    plot(Om,amplitude,'-','DisplayName',strrep(run_names{irun},'_','-'));
    figure(2)
    semilogy(Om,max(half_norm,1e-16),'-','DisplayName',strrep(run_names{irun},'_','-'));
end
for ifig=1:2
    figure(ifig)
    xlabel('$\omega_1$','Interpreter','latex');
    if ifig==1
        ylabel('$\|q\|_{\infty}$','Interpreter','latex');
    else
        set(gca,'YScale','log');
        ylabel('Half-integer coefficient norm');
    end
    grid on;
    legend('show','Location','best');
    hold off
    plot_setting_paper()
    legend boxon;
end
