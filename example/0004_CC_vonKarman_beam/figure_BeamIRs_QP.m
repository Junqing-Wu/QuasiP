% FIGURE_BEAMIRS_QP  Compare FSE shooting and mVCF at the paper's samples.
bd = BeamIRs_SHT_data();
[~,data] = QuasiP_read_solution('','beamIRs_po_stab',1);
dof = data.amp_opts.outdof(1);
om = coco_bd_col(bd,'Om');
amp = coco_bd_col(bd,'Max_x');
wave = coco_bd_col(bd,'MF_f',false);
coef = coco_bd_col(bd,'MF_x',false);
qpContent = cellfun(@(f,x) norm(x(f(:,4)~=0,dof)),wave,coef);
isQP = qpContent>1e-4;
assert(all(isQP),'The SHT branch includes points that collapsed to the periodic orbit.');

types = coco_bd_col(bd,'TYPE',false);
comparison = zeros(3,5);
for ii=1:3
    row = find(strcmp(types,sprintf('SHT%d',ii)));
    assert(isscalar(row),'The SHT sample labels are incomplete.');
    mvcf = coco_bd_read(sprintf('beamIRs_qp_mvcf_%d',ii));
    omM = coco_bd_col(mvcf,'Om');
    ampM = coco_bd_col(mvcf,'Max_x');
    fM = coco_bd_col(mvcf,'MF_f');
    xM = coco_bd_col(mvcf,'MF_x');
    assert(size(omM,1)==1 && abs(omM(1)-om(row,1))<1e-7, ...
        'The mVCF correction did not return the requested frequency.');
    assert(norm(xM(fM(:,4)~=0,dof))>1e-4,'The mVCF correction returned a periodic orbit.');
    comparison(ii,:) = [om(row,1),amp(row,dof),ampM(dof),om(row,2),omM(2)];
end
disp(array2table(comparison,'VariableNames', ...
    {'Omega1','Amp_SHT','Amp_mVCF','Omega2_SHT','Omega2_mVCF'}));

figure(1)
clf
plot(om(:,1),amp(:,dof),'-','Color',[0.85,0.35,0],'LineWidth',1.5)
hold on
plot(comparison(:,1),comparison(:,3),'mx','MarkerSize',8,'LineWidth',1.5)
xlabel('$\omega_1$ [rad/s]','Interpreter','latex');
ylabel('$\|w(l/4)\|_\infty$ [mm]','Interpreter','latex');
legend('FSE shooting','mVCF','Location','best');
grid on
plot_setting_paper()
hold off

figure(2)
clf
plot(om(:,1),om(:,2),'-','Color',[0.85,0.35,0],'LineWidth',1.5)
hold on
plot(comparison(:,1),comparison(:,5),'mx','MarkerSize',8,'LineWidth',1.5)
xlabel('$\omega_1$ [rad/s]','Interpreter','latex');
ylabel('$\omega_2$ [rad/s]','Interpreter','latex');
legend('FSE shooting','mVCF','Location','best');
grid on
plot_setting_paper()
hold off
