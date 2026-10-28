bd1 = coco_bd_read('qp511a_sys1');
max_x_a_1    = coco_bd_col(bd1, 'Max_x');
Om_model_a_1 = coco_bd_col(bd1, 'Om');
gamma_a_1    = coco_bd_col(bd1, 'gamma');
lab    = coco_bd_col(bd1,'LAB');
lab_vec = zeros(length(lab),1);
for ii = 1:length(lab)
    if ~isempty(lab{ii})
        lab_vec(ii) = lab{ii};
    else
        lab_vec(ii) = 0;
    end
end
labOm13  = coco_bd_labs(bd1, 'Om1=3');
idx_Om13 = ismember(lab_vec, labOm13);
idx_Om13_b = ismember(lab_vec, labOm13(4:5));



bd2 = coco_bd_read('qp511a_sys2');
max_x_a_2    = coco_bd_col(bd2, 'Max_x');
Om_model_a_2 = coco_bd_col(bd2, 'Om');
gamma_a_2    = coco_bd_col(bd2, 'gamma');

bd3 = coco_bd_read('qp511a_sys3');
max_x_a_3    = coco_bd_col(bd3, 'Max_x');
Om_model_a_3 = coco_bd_col(bd3, 'Om');
gamma_a_3    = coco_bd_col(bd3, 'gamma');


bd1 = coco_bd_read('qp511b_sys1');
max_x_b_1    = coco_bd_col(bd1, 'Max_x');
Om_model_b_1 = coco_bd_col(bd1, 'Om');
gamma_b_1  = coco_bd_col(bd1, 'gamma');

bd2 = coco_bd_read('qp511b_sys2');
max_x_b_2    = coco_bd_col(bd2, 'Max_x');
Om_model_b_2 = coco_bd_col(bd2, 'Om');
gamma_b_2  = coco_bd_col(bd2, 'gamma');

bd3 = coco_bd_read('qp511b_sys3');
max_x_b_3    = coco_bd_col(bd3, 'Max_x');
Om_model_b_3 = coco_bd_col(bd3, 'Om');
gamma_b_3  = coco_bd_col(bd3, 'gamma');



figure(1)
c_gray = [0.40, 0.40, 0.40];
c_orange = [0.85, 0.325, 0.098];
c_blue = [0.00, 0.447, 0.741];


plot(Om_model_a_1(:,1), max_x_a_1,'-','Color',c_gray,'LineWidth',1.2);
hold on
plot(Om_model_a_2(:,1), max_x_a_2(:,1),'--','Color',c_orange,'LineWidth',1.2);

[omSample_a,x3Sample_a,gammaSample_a,~] = ...
    sample_curve_points(Om_model_a_3(:,1),max_x_a_3(:,1),gamma_a_3,60,0.2);
plot(omSample_a, x3Sample_a,...
    'o','Color',c_blue,'MarkerSize',8);

plot(Om_model_a_1(idx_Om13,1), max_x_a_1(idx_Om13),'r^','MarkerSize',10);

xlabel('\omega_1');
ylabel('$\Vert x \Vert_{\infty}$','Interpreter','latex');
grid on;
legend('Sys. (1)','Sys. (2)','Sys. (3)','$\omega_1=3$','Interpreter','latex','Location','best')
xlim([1.8 5])
zlim([0 9])
plot_setting_paper()


figure(2)
c_gray = [0.40, 0.40, 0.40];
c_orange = [0.85, 0.325, 0.098];
c_blue = [0.00, 0.447, 0.741];

plot(gamma_b_1, max_x_b_1,'-','Color',c_gray,'LineWidth',1.2);
hold on
plot(gamma_b_2, max_x_b_2(:,1),'--','Color',c_orange,'LineWidth',1.2);

[gammaSample_b,x3Sample_b,omSample_b,~] = ...
    sample_curve_points(gamma_b_3,max_x_b_3(:,1),Om_model_b_3(:,1),35,0.2);
plot(gammaSample_b, x3Sample_b,...
    'o','Color',c_blue,'MarkerSize',8);
plot(gamma_a_1(idx_Om13_b), max_x_a_1(idx_Om13_b),'r^','MarkerSize',10);


xlabel('\gamma');
ylabel('$\Vert x \Vert_{\infty}$','Interpreter','latex');
grid on;
legend('Sys. (1)','Sys. (2)','Sys. (3)','$\omega_1=3$','Interpreter','latex','Location','best')
% xlim([1.8 5])
% zlim([0 9])
plot_setting_paper()