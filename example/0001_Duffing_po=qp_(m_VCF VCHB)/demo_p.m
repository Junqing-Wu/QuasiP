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
rho = 1;
Om1 = 4;
OM0 = rho*Om1;
pname = {'zeta' 'gamma' 'f1'};
p0 = [0.1;0.2;5];
%----------------------------------------------------------
% initial condition
Fex1.location{1} = 1;
Fex1.amplitude{1} = -1i*p0(3);
Fex1.location{2} = 50;
Fex1.amplitude{2} = -1i*p0(3);

[~,~,A,B,C] = PDuffing1_L(p0);
for ii = 1:2
    IP.amplitude{ii} = ( -(OM0)^2*A + 1i*OM0*B + C )\ ...
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
ode_opts    = struct('n',1, 'd',1, 'e',1, 'isfwithu',false);
amp_opts    = struct('isamp',true,'isuamp',true,'outdof',1,'Nperiod',1,'Nsampling',2^9);
stable_opts = struct('issta',true,'S',2^6,'L',2^8,'N_Ly',500);
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts,'stable_opts',stable_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',    'mVCF'           , ...     
    'type',  { {'HB'} }, ...   
    'k',       { {35} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^7} }, ...
    'S_wave',  { {2^8} });
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 

funcs = {@PDuffing1_F, @PDuffing1_DFDX, @PDuffing1_DFDP, @PDuffing1_L};
prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, X_MF2, OM0, pname, p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
bd1 = coco(prob1, 'qp11', [], 1, {'k_OM1' 'k_OM2' 'zeta' 'gamma' 'f1'}, [1.5 4]);
close(11111)

%  coco_bd_read('qp11');
max_x1    = coco_bd_col(bd1, 'Max_x');
Om_mode1  = coco_bd_col(bd1, 'Om');
Stab1      = coco_bd_col(bd1, 'Stab');

TX      = coco_bd_col(bd1, 'HT_x');
TU      = coco_bd_col(bd1, 'HT_u'); 
T = coco_bd_col(bd1, 't'); 


A1 = find(~Stab1);

figure(2)
plot(Om_mode1, max_x1,'k-');
hold on
plot(Om_mode1(A1), max_x1(A1), 'ko');
xlabel('\omega_1');
ylabel('max(x)');
grid on;
dd = 1;
figure(3)
plot(TX(:,dd),TU(:,dd))
hold on


