% ddx + 2*zeta*dx + x + gama*x^3 = f1*sin(Om1*t) + f1*sin(Om2*t)
% A(p)*\ddot{x} + D*\dot{x} + K*x = f(x,\dot{x},p,theta) 
%---------------------------------------------------------- 
% duffing1: A = 1; B = 2*zeta; C = 1; 
%           f = -gama*x^3 + f1*sin(theta1) + f1*sin(theta2)
%----------------------------------------------------------
% duffing2: A = []; B = [1 0;0 1]; C = []; 
%           f = [dx
%                -2*zeta*dx - x -gama*x^3 + f1*sin(theta1) + f1*sin(theta2)]
%----------------------------------------------------------
% duffing3: A = []; B = [2*zeta 1;1 0]; C = [1 0;0 -1]; 
%           f = [-gama*x^3 + f1*sin(theta1) + f1*sin(theta2)
%                0]
%----------------------------------------------------------
clearvars;
close all;
clc;
%% ====================== sys 1 ======================
bd1 = coco_bd_read('qp511a_sys1');
lab  = coco_bd_labs(bd1, 'Om1=3');
lab  = min(lab);
prob = coco_prob();
prob = ode_QuasiP2QuasiP(prob, '', 'qp511b_1', lab);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob, 'qp511b_sys1', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)
%% ====================== sys 2 ======================
bd1 = coco_bd_read('qp511a_sys2');
lab  = coco_bd_labs(bd1, 'Om1=3');
lab  = min(lab);
prob = coco_prob();
prob = ode_QuasiP2QuasiP(prob, '', 'qp511b_2', lab);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob, 'qp511b_sys2', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)
%% ====================== sys 3 ======================
bd1 = coco_bd_read('qp511a_sys3');
lab  = coco_bd_labs(bd1, 'Om1=3');
lab  = min(lab);
prob = coco_prob();
prob = ode_QuasiP2QuasiP(prob, '', 'qp511b_3', lab);
[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob = coco_add_func(prob, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
coco(prob, 'qp511b_sys3', [], 1, {'gamma' 'k_OM1' 'k_OM2' 'zeta' 'f1' 'f2'}, [0.1 0.6]);
close(11111)

function [data, y, J] = fc_2(prob, data, u)
Om1 = u(1);
Om2 = u(2);
y = Om1/2^0.5  - Om2;
J = [1/2^0.5 , -1];
end






