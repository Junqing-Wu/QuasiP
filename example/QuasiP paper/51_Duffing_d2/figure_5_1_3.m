bd1 = coco_bd_read('qp513_stability');
max_x1    = coco_bd_col(bd1, 'Max_x');
Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));

figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'b--','LineWidth',1.5}; % unstable
ST2{2} = {'b-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),max_x1,Stab1(:,1),ST2,'Unstable','Stable')
xlabel('\omega_1');
ylabel('$\Vert x \Vert_{\infty}$','Interpreter','latex');
grid on;
xlim([1.8 5])
zlim([0 9])
plot_setting_paper()



figure(2)
c_orange = [0.85, 0.325, 0.098];
c_blue = [0.00, 0.447, 0.741];

plot(332:-1:1,Ly1(:,1),'-','Color',c_orange,'LineWidth',1.2);
hold on
plot(332:-1:1,Ly1(:,3),'--','Color',c_blue,'LineWidth',1.2);
plot(332:-1:1,Ly1(:,2),'-','Color',c_orange,'LineWidth',1.2);
plot(332:-1:1,Ly1(:,4),'--','Color',c_blue,'LineWidth',1.2);

xlabel('Continuation number');
ylabel('$\sigma$','Interpreter','latex');
legend('Stab. 1','Stab. 2','Location','best')
grid on;
xlim([1 332])
plot_setting_paper()


