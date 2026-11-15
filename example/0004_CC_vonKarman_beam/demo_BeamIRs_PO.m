% DEMO_BEAMIRS_PO  Full-order periodic response of the clamped-pinned beam.
% Reference: "Parallelized computation of quasi-periodic solutions for finite
% element problems: A Fourier series expansion-based shooting method".
% Section 4.2, Fig. 15. Run COCO startup and QuasiP install first.
clearvars;
close all;
clc;

%% ====================== Beam model ======================
nElements = 32;
[M,C,K,fnl,E,output] = build_vonKarmanBeamIRs(nElements);
model = struct('M',M,'C',C,'K',K,'fnl',fnl,'E',full(E));
n = size(M,1);
outdof = [find(output(1,:)),find(output(2,:))];
assert(nElements==4*floor(nElements/4), ...
    'Use a multiple of four for the nodal quarter-span and midspan outputs.');

%% ====================== Initial periodic solution ======================
pname = {'force_scale'};
p0 = 1;
OM0 = 33;
% Linear complex response to cos(theta), using real(Xc*exp(i*theta)).
Xk_MF.k_wave = 1;
Xk_MF.Xc = ((K-OM0^2*M+1i*OM0*C)\(p0*E)).';

%% ====================== QuasiP continuation ======================
prob = coco_prob();
ode_opts = struct('n',n,'d',1,'e',1,'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',false,'isht',true, ...
    'outdof',outdof,'Nperiod',1,'Nsampling',256);
sta_opts = struct('issta',true,'L1',256,'Omt',1,'store_all_psi',false);
PHI_INITIAL = struct('param','MP','method','mVCF','isvar',true, ...
    'type',{{'HB'}},'k',{{7}},'P',{{1}},'S',{{32}},'S_wave',{{256}});
prob = coco_set(prob,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts, ...
    'sta_opts',sta_opts,'PHI_INITIAL',PHI_INITIAL);
prob = coco_set(prob,'cont','PtMX',[0,500],'NSV',5, ...
    'h0',0.025,'h_min',1e-5,'h_max',0.30,'MaxRes',50,'TOL',1e-10);
prob = coco_set(prob,'corr','ItMX',20,'TOL',1e-7);

funcs = {@(x,p,theta) BeamIRs_F(x,p,theta,model), ...
         @(x,p,theta) BeamIRs_DFDX(x,p,theta,model), ...
         @(x,p,theta) BeamIRs_DFDP(x,p,theta,model), ...
         @(p) BeamIRs_L(p,model), ...
         @(x,p,theta) BeamIRs_DFDTH(x,p,theta,model)};
bd = coco(ode_isol2QuasiP(prob,'',funcs{:},Xk_MF,OM0,pname,p0), ...
    'beamIRs_po_stab',[],1,{'k_OM1','force_scale','QuasiP.test.TR'},[33,34.5]);
