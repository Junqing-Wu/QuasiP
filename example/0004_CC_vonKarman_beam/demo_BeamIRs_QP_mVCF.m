% DEMO_BEAMIRS_QP_MVCF  Correct three SHT solutions by full-order mVCF.
% Reference: "Parallelized computation of quasi-periodic solutions for finite
% element problems: A Fourier series expansion-based shooting method".
% Section 4.2, Fig. 16. Run demo_BeamIRs_QP_SHT first.
clearvars;
close all;
clc;

%% ====================== Read the SHT sample labels ======================
[bd,runs] = BeamIRs_SHT_data();
saved = coco_bd_col(bd,'LAB',false);
types = coco_bd_col(bd,'TYPE',false);
targets = [33.859,33.931,34.034];
k1 = [11,13,11];

%% ====================== Fixed-frequency mVCF corrections ======================
for ii = 1:numel(targets)
    row = find(strcmp(types,sprintf('SHT%d',ii)));
    assert(isscalar(row),'The SHT branch must reach all sample events.');
    lab = saved{row};
    point = bd([1,row+1],:);
    [~,data] = QuasiP_read_solution('',runs{row},lab);
    OM0 = coco_bd_col(point,'Om');
    p0 = coco_bd_col(point,'p');
    assert(abs(OM0(1)-targets(ii))<1e-7,'SHT sample frequency mismatch.');
    MF_f = coco_bd_col(point,'MF_f');
    Xk_MF.k_wave = MF_f(:,3:4);
    Xk_MF.Xc = coco_bd_col(point,'MF_x');
    assert(size(Xk_MF.Xc,2)==data.ode_opts.n && ...
        norm(Xk_MF.Xc(Xk_MF.k_wave(:,2)~=0,:),'fro')>1e-3, ...
        'The SHT seed must carry nonzero second-frequency content in the full model.');
    % Project the SHT spectrum onto the retained HB modes before sampling.
    keep = abs(Xk_MF.k_wave(:,1))<=k1(ii) & abs(Xk_MF.k_wave(:,2))<=7;
    Xk_MF.k_wave = Xk_MF.k_wave(keep,:);
    Xk_MF.Xc = Xk_MF.Xc(keep,:);

    prob = coco_prob();
    amp_opts = struct('isamp',true,'isuamp',false,'isht',false, ...
        'outdof',data.amp_opts.outdof,'Nperiod',0,'Nsampling',256);
    sta_opts = struct('issta',false);
    PHI_INITIAL = struct('param','MP','method','mVCF','type',{{'HB','HB'}}, ...
        'isvar',true,'sparse_jacobian',true,'k',{{k1(ii),7}}, ...
        'P',{{1,1}},'S',{{32,32}},'S_wave',{{128,32}});
    prob = coco_set(prob,'QuasiP','ode_opts',data.ode_opts,'amp_opts',amp_opts, ...
        'sta_opts',sta_opts,'PHI_INITIAL',PHI_INITIAL);
    prob = coco_set(prob,'corr','ItMX',15,'SubItMX',1,'TOL',1e-7);
    funcs = {data.f,data.fx,data.fp,data.l,data.fth};
    prob = ode_isol2QuasiP(prob,'',funcs{:},Xk_MF,OM0,data.pnames_np,p0);
    run = sprintf('beamIRs_qp_mvcf_%d',ii);
    corrected = coco(prob,run,[],0,{'u_om2','k_OM1','force_scale'});
    assert(isscalar(coco_bd_labs(corrected,'EP')), ...
        'The mVCF corrector failed at omega1=%.3f.',targets(ii));
end

