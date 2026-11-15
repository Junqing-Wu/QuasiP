%% ====================== Plot tor11 ======================
bd_tor11   = coco_bd_read('qp515_tor11');
HT_x_tor11 = coco_bd_col(bd_tor11,'HT_x');
HT_u_tor11 = coco_bd_col(bd_tor11,'HT_u');
T_x_tor11  = coco_bd_col(bd_tor11,'T_x');
T_u_tor11  = coco_bd_col(bd_tor11,'T_u');
Om_tor11   = coco_bd_col(bd_tor11, 'Om');

bd_tor11   = coco_bd_read('qp515_tor11_b');
T_x_tor11_b  = coco_bd_col(bd_tor11,'T_x');
T_u_tor11_b  = coco_bd_col(bd_tor11,'T_u');


figure(11)
plot(x_tor11(1:end/8,1),x_tor11(1:end/8,2), 'color', [0.5 0.5 0.5])
hold on
plot(HT_x_tor11,HT_u_tor11,'r.','MarkerSize',5)
plot(T_x_tor11(3:2^2:end),T_u_tor11(3:2^2:end),'g.','MarkerSize',4)
plot(T_x_tor11_b(2:2^2:end),T_u_tor11_b(2:2^2:end),'c.','MarkerSize',4)
legend('NI','QuasiP','$x(2k\pi+\pi,\tau_2)$','$x(\tau_1,2k\pi+\pi/2)$','Location','best')
xlim([-10,10])
ylim([-25,25])
xlabel('$x$','Interpreter','latex');
ylabel('$\dot{x}$','Interpreter','latex');
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)

MF_x_tor11 = coco_bd_col(bd_tor11,'MF_x');
MF_f_tor11 = coco_bd_col(bd_tor11,'MF_f');
figure(12)
plot(ff_tor11/(Om_tor11(1)),20*log10(abs(X_tor11(:,1))),'k-')
hold on
plot(MF_f_tor11(:,1),20*log10(abs(2*MF_x_tor11)),'ro','MarkerSize',8)
xlabel('$f$','Interpreter','latex')
ylabel('$\mathrm{dB}$','Interpreter','latex')
xlim([0,6])
ylim([-125,25])
legend('NI','QuasiP','Location','best')
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)

%% ====================== Plot tor22 ======================
bd_tor22   = coco_bd_read('qp515_tor22');
HT_x_tor22 = coco_bd_col(bd_tor22,'HT_x');
HT_u_tor22 = coco_bd_col(bd_tor22,'HT_u');
T_x_tor22  = coco_bd_col(bd_tor22,'T_x');
T_u_tor22  = coco_bd_col(bd_tor22,'T_u');
Om_tor22   = coco_bd_col(bd_tor22, 'Om');

bd_tor22   = coco_bd_read('qp515_tor22_b');
T_x_tor22_b  = coco_bd_col(bd_tor22,'T_x');
T_u_tor22_b  = coco_bd_col(bd_tor22,'T_u');

figure(21)
plot(x_tor22(1:end/8,1),x_tor22(1:end/8,2), 'color', [0.5 0.5 0.5])
hold on
plot(HT_x_tor22,HT_u_tor22,'r.','MarkerSize',5)
plot(T_x_tor22(3:2^2:end),T_u_tor22(3:2^2:end),'g.','MarkerSize',4)
plot(T_x_tor22_b(2:2^2:end),T_u_tor22_b(2:2^2:end),'c.','MarkerSize',4)
legend('NI','QuasiP','$x(2k\pi+\pi,\tau_2)$','$x(\tau_1,2k\pi+\pi/2)$','Location','best')
xlim([-10,10])
ylim([-25,25])
xlabel('$x$','Interpreter','latex');
ylabel('$\dot{x}$','Interpreter','latex');
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)

MF_x_tor22 = coco_bd_col(bd_tor22,'MF_x');
MF_f_tor22 = coco_bd_col(bd_tor22,'MF_f');
figure(22)
plot(ff_tor22/(Om_tor22(1)),20*log10(abs(X_tor22(:,1))),'k-')
hold on
plot(MF_f_tor22(:,1),20*log10(abs(2*MF_x_tor22)),'ro','MarkerSize',8)
xlabel('$f$','Interpreter','latex')
ylabel('$\mathrm{dB}$','Interpreter','latex')
legend('NI','QuasiP','Location','best')
xlim([0,6])
ylim([-125,25])
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)

%% ====================== Plot tor44 ======================
bd_tor44   = coco_bd_read('qp515_tor44');
HT_x_tor44 = coco_bd_col(bd_tor44,'HT_x');
HT_u_tor44 = coco_bd_col(bd_tor44,'HT_u');
T_x_tor44  = coco_bd_col(bd_tor44,'T_x');
T_u_tor44  = coco_bd_col(bd_tor44,'T_u');
Om_tor44   = coco_bd_col(bd_tor44, 'Om');

bd_tor44   = coco_bd_read('qp515_tor44_b');
T_x_tor44_b  = coco_bd_col(bd_tor44,'T_x');
T_u_tor44_b  = coco_bd_col(bd_tor44,'T_u');

figure(41)
plot(x_tor44(1:end/8,1),x_tor44(1:end/8,2), 'color', [0.5 0.5 0.5])
hold on
plot(HT_x_tor44,HT_u_tor44,'r.','MarkerSize',5)
plot(T_x_tor44(3:2^2:end),T_u_tor44(3:2^2:end),'g.','MarkerSize',4)
plot(T_x_tor44_b(2:2^2:end),T_u_tor44_b(2:2^2:end),'c.','MarkerSize',4)
legend('NI','QuasiP','$x(2k\pi+\pi,\tau_2)$','$x(\tau_1,2k\pi+\pi/2)$','Location','best')
xlim([-10,10])
ylim([-25,25])
xlabel('$x$','Interpreter','latex');
ylabel('$\dot{x}$','Interpreter','latex');
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)

MF_x_tor44 = coco_bd_col(bd_tor44,'MF_x');
MF_f_tor44 = coco_bd_col(bd_tor44,'MF_f');
figure(42)
plot(ff_tor44/(Om_tor44(1)),20*log10(abs(X_tor44(:,1))),'k-')
hold on
plot(MF_f_tor44(:,1),20*log10(abs(2*MF_x_tor44)),'ro','MarkerSize',8)
xlabel('$f$','Interpreter','latex')
ylabel('$\mathrm{dB}$','Interpreter','latex')
xlim([0,6])
ylim([-125,25])
legend('NI','QuasiP','Location','best')
hold off
plot_setting_paper()
set(findall(gcf,'type','line'),'linewidth',1)
