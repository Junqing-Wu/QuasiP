bd1 = coco_bd_read('qp512_meth1');
t1    = coco_bd_col(bd1, 't');
x1  = coco_bd_col(bd1, 'T_x');
u1  = coco_bd_col(bd1, 'T_u');

bd2 = coco_bd_read('qp512_meth2');
t2    = coco_bd_col(bd2, 't');
x2  = coco_bd_col(bd2, 'T_x');
u2  = coco_bd_col(bd2, 'T_u');

bd3 = coco_bd_read('qp512_meth3');
t3    = coco_bd_col(bd3, 't');
x3  = coco_bd_col(bd3, 'T_x');
u3  = coco_bd_col(bd3, 'T_u');

bd4 = coco_bd_read('qp512_meth4');
t4    = coco_bd_col(bd4, 't');
x4  = coco_bd_col(bd4, 'T_x');
u4  = coco_bd_col(bd4, 'T_u');

bd5 = coco_bd_read('qp512_meth5');
t5    = coco_bd_col(bd5, 't');
x5  = coco_bd_col(bd5, 'T_x');
u5  = coco_bd_col(bd5, 'T_u');


c1 = [0.40, 0.40, 0.40];
c2 = [0.85, 0.325, 0.098];
c3 = [0.00, 0.447, 0.741];
c4 = [0.494, 0.184, 0.556];
c5 = [0.466, 0.674, 0.188];

figure(1)
plot(t1,x1,'-','Color',c1,'LineWidth',1.2);
hold on
plot(t2,x2,'--','Color',c2,'LineWidth',1.2);
plot(t3(1:8:end),x3(1:8:end),'s','Color',c3,'MarkerSize',6);
plot(t4(1:8:end),x4(1:8:end),'^','Color',c4,'MarkerSize',6);
plot(t5(1:8:end),x5(1:8:end),'o','Color',c5,'MarkerSize',6);
hold off
xlabel('$t$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
legend('Meth. 1','Meth. 2','Meth. 3','Meth. 4','Meth. 5');
plot_setting_paper()

figure(2)
plot(t1,u1,'-','Color',c1,'LineWidth',1.2);
hold on
plot(t2,u2,'--','Color',c2,'LineWidth',1.2);
plot(t3(1:8:end),u3(1:8:end),'s','Color',c3,'MarkerSize',6);
plot(t4(1:8:end),u4(1:8:end),'^','Color',c4,'MarkerSize',6);
plot(t5(1:8:end),u5(1:8:end),'o','Color',c5,'MarkerSize',6);
hold off
xlabel('$t$','Interpreter','latex');
ylabel('$\dot{x}$','Interpreter','latex');
legend('Meth. 1','Meth. 2','Meth. 3','Meth. 4','Meth. 5');
plot_setting_paper()

