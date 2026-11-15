% DEMO_5_1_1B  Restart the three formulations in cubic stiffness gamma.
% Article: Section 5.1, Fig. 6(b).
% Requires qp511a_sys1/2/3 from demo_5_1_1a and saved Om1=3 events.
% Select the smallest available event label for each source run.
% Continue gamma in [0.1,0.6]; keep the prescribed excitation frequencies fixed.
% Outputs: qp511b_sys1, qp511b_sys2, qp511b_sys3.
% Next: figure_5_1_1.
% Setup: execute COCO startup and QuasiP install; work in this example folder.
% COCO outputs normally go to data/<run_name>/. See README.md for the running order.
% --------------------------------------------------------------------------
% Two-frequency Duffing oscillator (QuasiP article, Sec. 5.1).
% Physical displacement q(t) satisfies
%   qddot + 2*zeta*qdot + q + gamma*q^3
%       = f1*sin(Om1*t) + f2*sin(Om2*t).
% Here theta1 = Om1*t, theta2 = Om2*t, and p = [zeta; gamma; f1; f2].
% The implemented linear/nonlinear split is
%   A(p)*xddot + B(p)*xdot + C(p)*x = f(x,p,theta).
% -------------------------------------------------------------------------
% Sys. (1), Duffing1: x = q (scalar second-order formulation).
%   A = 1; B = 2*zeta; C = 1;
%   f = -gamma*x^3 + f1*sin(theta1) + f2*sin(theta2).
% -------------------------------------------------------------------------
% Sys. (2), Duffing2: x = [x1; x2] = [q; qdot].
%   A is absent; B = eye(2); C is absent (no L callback is required).
%   xdot = [x2;
%           -2*zeta*x2 - x1 - gamma*x1^3 ...
%           + f1*sin(theta1) + f2*sin(theta2)].
% -------------------------------------------------------------------------
% Sys. (3), Duffing3: x = [x1; x2] = [q; qdot].
%   A is absent; B = [2*zeta, 1; 1, 0]; C = [1, 0; 0, -1];
%   B*xdot + C*x = [-gamma*x1^3 + f1*sin(theta1) + f2*sin(theta2); 0].
%   The second row enforces xdot1 = x2.
% -------------------------------------------------------------------------
clearvars;
close all;
clc;
%% ====================== sys 1 ======================
% Read prerequisite saved branch data; the corresponding producer demo must run first.
bd1 = coco_bd_read('qp511a_sys1');
lab  = coco_bd_labs(bd1, 'Om1=3');
lab  = min(lab);
prob = coco_prob();
% Restore the selected saved solution as the starting point of a new continuation problem.
prob = ode_QuasiP2QuasiP(prob, '', 'qp511a_sys1', lab);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob, 'qp511b_sys1', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)
%% ====================== sys 2 ======================
% Read prerequisite saved branch data; the corresponding producer demo must run first.
bd1 = coco_bd_read('qp511a_sys2');
lab  = coco_bd_labs(bd1, 'Om1=3');
lab  = min(lab);
prob = coco_prob();
% Restore the selected saved solution as the starting point of a new continuation problem.
prob = ode_QuasiP2QuasiP(prob, '', 'qp511a_sys2', lab);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob, 'qp511b_sys2', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)
%% ====================== sys 3 ======================
% Read prerequisite saved branch data; the corresponding producer demo must run first.
bd1 = coco_bd_read('qp511a_sys3');
lab  = coco_bd_labs(bd1, 'Om1=3');
lab  = min(lab);
prob = coco_prob();
% Restore the selected saved solution as the starting point of a new continuation problem.
prob = ode_QuasiP2QuasiP(prob, '', 'qp511a_sys3', lab);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
% Track the one-dimensional solution branch (COCO dimension 1).
% The requested parameter bounds do not guarantee actual coverage; inspect the run log.
coco(prob, 'qp511b_sys3', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)

% Zero constraint preserving the prescribed excitation-frequency ratios.
% Return its analytic Jacobian with respect to the supplied frequency vector.
function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end






