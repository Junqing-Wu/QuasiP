function plot_setting_paper()
%PLOT_SETTING_PAPER  Apply consistent figure settings for QuasiP examples.
%
% Set the fonts, line widths, axes appearance, renderer, and other graphics
% defaults used by the paper and example plotting scripts.
%
% On input:
%   None.
%
% On output:
%   None. MATLAB graphics defaults are updated.

    set(gca,'linewidth',1.5);
    set(gca,'box','on','linewidth',1.5);
    set(gca,'FontName','Times New Roman','fontsize',30);
    legend('FontSize',25,'Interpreter','latex')
    set(findall(gcf,'type','line'),'linewidth',2)
%     axis tight
    set(gcf,'unit','normalized','position',[0.1,0.1,0.32,0.48]);
    legend boxoff;
end