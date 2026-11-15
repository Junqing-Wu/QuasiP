% DEMO_5_1_4C  Continue Tor44 from a saved representative solution.
% Article: Section 5.1, Fig. 9; P=(4,4), still two base frequencies.
% Requires qp515_tor44 from demo_5_1_5; initial Om1=2.71.
% Reconstruct the initial Fourier guess from the saved MF_f and MF_x columns.
% Outputs: qp514_tor44_a and qp514_tor44_b.
% Next: figure_5_1_4 after Tor11 and all Tor22 branches are available.
% Setup: execute COCO startup and QuasiP install; work in this example folder.
% COCO outputs normally go to data/<run_name>/. See README.md for the running order.
% Model: qddot + 2*zeta*qdot + q + gamma*q^3
%        = f1*sin(Om1*t) + f2*sin(Om2*t); theta_i=Om_i*t.
% --------------------------------------------------------------------------
%% ====================== tor44 ======================
% Define excitation frequencies and model parameters for this experiment.
rho = [1;1/2^0.5];
Om1 = 2.71;%2.71----2.75----2.85
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% Prepare the initial torus guess (integration or saved Fourier data below).
% Direct integration is used only when present; otherwise read a saved torus.
% Read prerequisite saved branch data; the corresponding producer demo must run first.
bd_tor44   = coco_bd_read('qp515_tor44');
MF_f       = coco_bd_col(bd_tor44,'MF_f');
MF_x       = coco_bd_col(bd_tor44,'MF_x');
Xk_MF2_tor44.k_wave = MF_f(:,3:4);
Xk_MF2_tor44.Xc = MF_x;
%----------------------------------------------------------
%% ====================== tor44 ======================
prob1 = coco_prob();
% ODE settings: n is the number of model coordinates; d/e count response/excitation frequencies.
% isfwithu indicates whether the nonlinear callback explicitly depends on velocity.
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
% Amplitude settings control selected outputs, velocity reconstruction, and time-domain sampling.
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.5);
% Stability settings: S samples transverse phases, L1 resolves each trajectory,
% N_Ly sets QR iteration count, and Omt selects integration direction(s).
% Near-zero exponents require convergence checks; the stored binary flag alone is insufficient.
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
% Discretization: MP/mVCF approximates the full torus; MT/SHT integrates trajectory families.
% type selects HB/FD/CO or HT/MF; k sets basis order/stencil; P sets period factors.
% S is the solution sampling grid; S_wave is the denser reconstruction grid.
% OmOd sets MP frequency ordering; Omt and M select SHT integration direction and segment count.
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {4,4} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^7,2^7} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
% Model callbacks supply nonlinear force, state/parameter derivatives, and the linear part.
% Forced shooting additionally requires derivatives with respect to excitation phases.
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob1, 'qp514_tor44_a', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [2.3236 2.71]);


prob1 = coco_prob();
ode_opts = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',5000,'Nsampling',2^2);
prob1 = coco_set(prob1,'cont','PtMX',1000);
prob1 = coco_set(prob1,'cont','NSV',1);
prob1 = coco_set(prob1,'corr','ItMX',50);
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.25);
sta_opts = struct('issta',true,'S',2^5,'L1',2^8,'N_Ly',2000,'Omt',[1,2]);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'sta_opts',sta_opts);
PHI_INITIAL = struct( 'param', 'MP', ...
'method',   'mVCF'           , ...
'type',    { {'HB','HB'} }, ...
'k',       { {7,7} }, ...
'P',       { {4,4} }, ...
'S',       { {2^5,2^5} }, ...
'S_wave',  { {2^7,2^7} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL);
funcs = {@Duffing1_F, @Duffing1_DFDX, @Duffing1_DFDP, @Duffing1_L};
% Convert the Fourier initial guess to method-specific coefficients and register the COCO problem.
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2_tor44, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob1, 'qp514_tor44_b', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [2.71 2.72]);


% Zero constraint preserving the prescribed excitation-frequency ratios.
% Return its analytic Jacobian with respect to the supplied frequency vector.
function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end
