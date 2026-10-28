function plot_setting_paper()
    grid on
    set(gca,'linewidth',1.5);
    set(gca,'box','on','linewidth',1.5);
    set(gca,'FontName','Times New Roman','fontsize',30);
    legend('FontSize',25,'Interpreter','latex')
    set(findall(gcf,'type','line'),'linewidth',2)
%     axis tight
    set(gcf,'unit','normalized','position',[0.1,0.1,0.32,0.48]);
    legend boxoff;
end