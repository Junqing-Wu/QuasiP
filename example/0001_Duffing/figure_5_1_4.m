SS = 2^8;

PP = 1;
bd1 = coco_bd_read('qp514_tor11');
max_x1    = coco_bd_col(bd1, 'Max_x');
HT_x1     = coco_bd_col(bd1, 'HT_x');
HT_xa     = reshape(HT_x1,PP*SS,PP*SS,[]);
HT_xb     = pagetranspose(HT_xa);

HT_xa     = reshape(HT_xa(:,1,:),PP*SS,[]);
HT_xa     = transpose(HT_xa(129:SS:end,:));

HT_xb     = reshape(HT_xb(:,1,:),PP*SS,[]);
HT_xb     = transpose(HT_xb(65:SS:end,:));

Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));



figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'k:','LineWidth',1.5}; % unstable
ST2{2} = {'k-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xa,Stab1(:,1),ST2,'tor11 Unstable','tor11 Stable')
xlabel('\omega_1');
ylabel('$x(2k\pi+\pi,0)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'k:','LineWidth',1.5}; % unstable
ST2{2} = {'k-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xb,Stab1(:,1),ST2,'tor11 Unstable','tor11 Stable')
xlabel('\omega_1');
ylabel('$x(0,2k\pi+\pi/2)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()



PP = 2;
bd1 = coco_bd_read('qp514_tor22_b');
max_x1    = coco_bd_col(bd1, 'Max_x');
HT_x1     = coco_bd_col(bd1, 'HT_x');
HT_xa     = reshape(HT_x1,PP*SS,PP*SS,[]);
HT_xb     = pagetranspose(HT_xa);

HT_xa     = reshape(HT_xa(:,1,:),PP*SS,[]);
HT_xa     = transpose(HT_xa(129:SS:end,:));

HT_xb     = reshape(HT_xb(:,1,:),PP*SS,[]);
HT_xb     = transpose(HT_xb(65:SS:end,:));

Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));



figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,1),Stab1(:,1),ST2,'tor22 Unstable','tor22 Stable')
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,2),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(2k\pi+\pi,0)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,1),Stab1(:,1),ST2,'tor22 Unstable','tor22 Stable')
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,2),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(0,2k\pi+\pi/2)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


bd1 = coco_bd_read('qp514_tor22_bb');
max_x1    = coco_bd_col(bd1, 'Max_x');
HT_x1     = coco_bd_col(bd1, 'HT_x');
HT_xa     = reshape(HT_x1,PP*SS,PP*SS,[]);
HT_xb     = pagetranspose(HT_xa);

HT_xa     = reshape(HT_xa(:,1,:),PP*SS,[]);
HT_xa     = transpose(HT_xa(129:SS:end,:));

HT_xb     = reshape(HT_xb(:,1,:),PP*SS,[]);
HT_xb     = transpose(HT_xb(65:SS:end,:));

Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));



figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,1),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,2),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(2k\pi+\pi,0)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,1),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,2),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(0,2k\pi+\pi/2)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


PP = 2;
bd1 = coco_bd_read('qp514_tor22_a');
max_x1    = coco_bd_col(bd1, 'Max_x');
HT_x1     = coco_bd_col(bd1, 'HT_x');
HT_xa     = reshape(HT_x1,PP*SS,PP*SS,[]);
HT_xb     = pagetranspose(HT_xa);

HT_xa     = reshape(HT_xa(:,1,:),PP*SS,[]);
HT_xa     = transpose(HT_xa(129:SS:end,:));

HT_xb     = reshape(HT_xb(:,1,:),PP*SS,[]);
HT_xb     = transpose(HT_xb(65:SS:end,:));


Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));



figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,1),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,2),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(2k\pi+\pi,0)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'r:','LineWidth',1.5}; % unstable
ST2{2} = {'r-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,1),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,2),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(0,2k\pi+\pi/2)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()



PP = 4;
SS = 2^7;
bd1 = coco_bd_read('qp514_tor44_a');
max_x1    = coco_bd_col(bd1, 'Max_x');
HT_x1     = coco_bd_col(bd1, 'HT_x');
HT_xa     = reshape(HT_x1,PP*SS,PP*SS,[]);
HT_xb     = pagetranspose(HT_xa);

HT_xa     = reshape(HT_xa(:,1,:),PP*SS,[]);
HT_xa     = transpose(HT_xa(65:SS:end,:));

HT_xb     = reshape(HT_xb(:,1,:),PP*SS,[]);
HT_xb     = transpose(HT_xb(33:SS:end,:));

Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));



figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'b:','LineWidth',1.5}; % unstable
ST2{2} = {'b-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,1),Stab1(:,1),ST2,'tor44 Unstable','tor44 Stable')
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,2),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,3),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,4),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(2k\pi+\pi,0)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'b:','LineWidth',1.5}; % unstable
ST2{2} = {'b-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,1),Stab1(:,1),ST2,'tor44 Unstable','tor44 Stable')
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,2),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,3),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,4),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(0,2k\pi+\pi/2)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


bd1 = coco_bd_read('qp514_tor44_b');
max_x1    = coco_bd_col(bd1, 'Max_x');
HT_x1     = coco_bd_col(bd1, 'HT_x');
HT_xa     = reshape(HT_x1,PP*SS,PP*SS,[]);
HT_xb     = pagetranspose(HT_xa);

HT_xa     = reshape(HT_xa(:,1,:),PP*SS,[]);
HT_xa     = transpose(HT_xa(65:SS:end,:));

HT_xb     = reshape(HT_xb(:,1,:),PP*SS,[]);
HT_xb     = transpose(HT_xb(33:SS:end,:));

Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1     = coco_bd_col(bd1, 'Stab');
Ly1       = coco_bd_col(bd1, 'Ly');
A11 = find(~Stab1(:,1));



figure(1)
hold on
ST2 = cell(2,1);
ST2{1} = {'b:','LineWidth',1.5}; % unstable
ST2{2} = {'b-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,1),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,2),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,3),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xa(:,4),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(2k\pi+\pi,0)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()


figure(2)
hold on
ST2 = cell(2,1);
ST2{1} = {'b:','LineWidth',1.5}; % unstable
ST2{2} = {'b-','LineWidth',1.5};  % stable
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,1),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,2),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,3),Stab1(:,1),ST2)
plot_stab_lines4(Om_mode1(:,1),HT_xb(:,4),Stab1(:,1),ST2)
xlabel('\omega_1');
ylabel('$x(0,2k\pi+\pi/2)$','Interpreter','latex');
grid on;
% xlim([1.8 3])
% zlim([0 9])
plot_setting_paper()