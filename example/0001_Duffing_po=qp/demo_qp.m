% ddx + 2*zeta*dx + x + gama*x^3 = f1*sin(Om1*t) + f1*sin(Om2*t)
% A(p)*\ddot{x} + D*\dot{x} + K*x = f(x,\dot{x},p,theta) 
%---------------------------------------------------------- 
% duffing1: A = 1; B = 2*zeta; C = 1; 
%           f = -gama*x^3 + f1*sin(theta1) + f1*sin(RHO_2*theta1)
%----------------------------------------------------------
clearvars;
close all;
clc;
% system parameters
rho = [1;7];
Om1 = 4;
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1' 'f2'};
p0 = [0.1;0.2;5;5];
%----------------------------------------------------------
% initial condition
Fex1.location{1} = [1;0];
Fex1.amplitude{1} = -1i*p0(3);
Fex1.location{2} = [0;1];
Fex1.amplitude{2} = -1i*p0(4);

[~,~,A,B,C] = QPDuffing1_L(p0);
p = length(OM0);
for ii = 1:p
    IP.amplitude{ii} = ( -(OM0(ii))^2*A + 1i*OM0(ii)*B + C )\ ...
       Fex1.amplitude{ii} ;
    IP.location{ii} = Fex1.location{ii};
end
X_MF2 = [IP.location{1}',transpose(IP.amplitude{1})];
for ii = 2:numel(IP.location)
X_MF2 = [X_MF2;IP.location{ii}',transpose(IP.amplitude{ii})];
end
%----------------------------------------------------------

%% ====================== 任务 1 ======================
prob1 = coco_prob();
ode_opts    = struct('n',1, 'd',2, 'e',2, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',1,'Nsampling',2^12);
stable_opts = struct('issta',true,'S',2^5,'L',2^6,'N_Ly',500);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB'} }, ...   
    'k',       { {5,5} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^5,2^5} }, ...
    'S_wave',  { {2^7,2^7} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@QPDuffing1_F, @QPDuffing1_DFDX, @QPDuffing1_DFDP, @QPDuffing1_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
fc_funcs = {@fc_2};
prob1 = coco_add_func(prob1, 'KOM_zero', fc_funcs{:}, [], 'zero','uidx', uidx(data.kOM_idx) , 'F+dF');
bd1 =  coco(prob1, 'qp1', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1' 'f2'}, [1.5 4]);

% bd1 = coco_bd_read('qp1');
max_x1    = coco_bd_col(bd1, 'Max_x');
Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1      = coco_bd_col(bd1, 'Stab');
A1 = find(~Stab1);

figure(1)
plot(Om_mode1(:,1), max_x1,'b-');
hold on
plot(Om_mode1(A1,1), max_x1(A1), 'bo');
xlabel('\omega_1');
ylabel('max(x)');
grid on;

data = [Om_mode1(:,1), max_x1',Stab1'];
HT_x    = coco_bd_col(bd1, 'HT_x');
HT_u    = coco_bd_col(bd1, 'HT_u');
T_x    = coco_bd_col(bd1, 'T_x');
T_u    = coco_bd_col(bd1, 'T_u');
T    = coco_bd_col(bd1, 't');

dd = 1;
figure(2)
plot(HT_x(:,dd),HT_u(:,dd),'r.')
hold on
plot(T_x(:,dd),T_u(:,dd),'b-')




function [data, y, J] = fc_2(prob, data, u)
RHO_2 = 7;
Om1 = u(1);
Om2 = u(2);
y = RHO_2* Om1 - Om2;
J = [RHO_2 , -1];
end

