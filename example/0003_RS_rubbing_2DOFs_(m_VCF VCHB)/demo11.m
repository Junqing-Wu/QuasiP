clc;
clear;

bd22 = coco_bd_read('qp2');

SMF_x = coco_bd_col(bd22,'SX.MF_x');
SMF_f = coco_bd_col(bd22,'SX.MF_f');
SX_MF  = [SMF_f{1}(:,3:4),SMF_x{1}];

USMF_x = coco_bd_col(bd22,'USX.MF_x');
USMF_f = coco_bd_col(bd22,'USX.MF_f');
USX_MF  = [USMF_f{1}(:,3:4),USMF_x{1}];

SOM = coco_bd_col(bd22,'SX.Om');
USOM = coco_bd_col(bd22,'USX.Om');

SOM0    = SOM(1,:)';
USOM0   = USOM(1,:)';


SC  = coco_bd_col(bd22,'SX.C');
USC = coco_bd_col(bd22,'USX.C');
DC  = abs(norm(reshape(USC(:,:,1),[],1))-norm(reshape(SC(:,:,1),[],1)));

pname = {'mu','OM'};
p0 = [0.2;SOM0(1)];

prob2 = coco_prob();
ode_opts = struct('n',2, 'd',2, 'e',1, 'isfwithu',true);
amp_opts = struct('isamp',true,'isuamp',true,'outdof',[1 2],'Nperiod',50,'Nsampling',2^8);
prob2 = coco_set(prob2,'cont','PtMX',3000);
prob2 = coco_set(prob2,'cont','h_min',1e-3,'h_max',0.05);
prob2 = coco_set(prob2,'cont','NSV',1);
prob2 = coco_set(prob2,'QuasiP','ode_opts',ode_opts,'amp_opts',amp_opts); 
PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ...     
    'method',  { {'HB','HB'} }, ...   
    'k',       { {9,9} }, ...        
    'P',       { {1,1} }, ...    
    'S',       { {2^8,2^8} }, ...
    'S_wave',  { {2^8,2^8} });
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@RS1_F, @RS1_DFDX, @RS1_DFDP, @RS1_L};
prob2 = ode_isol2QuasiP(prob2, 'SX', funcs{:}, SX_MF, SOM0, pname , p0);
[SXdata, SXuidx] = coco_get_func_data(prob2, 'SX.QuasiP', 'data', 'uidx');
prob2 = coco_add_glue(prob2,'SX.OMp2',SXuidx(SXdata.kOM_idx),SXuidx(SXdata.p_idx(2)));

prob2 = ode_isol2QuasiP(prob2, 'USX', funcs{:}, USX_MF, USOM0, pname , p0);
[USXdata, USXuidx] = coco_get_func_data(prob2, 'USX.QuasiP', 'data', 'uidx');
prob2 = coco_add_glue(prob2,'USX.OMp2',USXuidx(USXdata.kOM_idx),USXuidx(USXdata.p_idx(2)));

% prob2 = coco_add_glue(prob2,'pname1',SXuidx(SXdata.p_idx(1)),USXuidx(USXdata.p_idx(1)));
% prob2 = coco_add_glue(prob2,'pname2',SXuidx(SXdata.p_idx(2)),USXuidx(USXdata.p_idx(2)));
prob2 = coco_add_glue(prob2,'KOM1',  SXuidx(SXdata.kOM_idx),USXuidx(USXdata.kOM_idx));

C_idx = [SXuidx(SXdata.C_idx);USXuidx(USXdata.C_idx)];

prob2 = coco_add_func(prob2,'my_dc',@my_DC,[],'inactive','DC','uidx',C_idx);

% bd2 = coco(prob2, 'qp2', [], 1, {'SX.k_OM1' 'Dom' 'SX.u_om2' 'USX.u_om2'  'USX.k_OM1' 'SX.OM' 'USX.OM' 'SX.mu' 'USX.mu'}, [0.3 0.45]);Dom2

coco(prob2, 'qp22', [], 1, {'DC' 'SX.u_om2' 'USX.u_om2' 'SX.k_OM1' 'USX.k_OM1' 'SX.OM' 'USX.OM' 'SX.mu' 'USX.mu' }, [0.04 0.1]);
close(11111)

bd22 = coco_bd_read('qp22');
bd1 = coco_bd_read('qp1');
OM1 = coco_bd_col(bd1,'k_OM1');
om2 = coco_bd_col(bd1,'u_om2');


SOM1 = coco_bd_col(bd22,'SX.k_OM1');
Som2 = coco_bd_col(bd22,'SX.u_om2');
USOM1 = coco_bd_col(bd22,'USX.k_OM1');
USom2 = coco_bd_col(bd22,'USX.u_om2');

max_x = coco_bd_col(bd1, 'Max_x');
Smax_x = coco_bd_col(bd22, 'SX.Max_x');
USmax_x = coco_bd_col(bd22, 'USX.Max_x');


figure(1)
plot(OM1,om2,'r--');
hold on
plot(SOM1,Som2,'bo');
plot(USOM1,USom2,'k^');
xlabel('k_{OM1}');
ylabel('\omega_2');
grid on;

figure(2)
plot(OM1,max_x(:,1),'r--');
hold on
plot(SOM1,Smax_x(:,1),'bo');
plot(USOM1,USmax_x(:,1),'k^');
xlabel('k_{OM1}');
ylabel('max(x)');
grid on;

Sht_x  = coco_bd_col(bd22, 'SX.HT_x');
USht_x = coco_bd_col(bd22, 'USX.HT_x');
Sht_u  = coco_bd_col(bd22, 'SX.HT_u');
USht_u = coco_bd_col(bd22, 'USX.HT_u');

iii = 1;
figure(3)
plot3(Sht_x(:,1,iii),Sht_x(:,2,iii),Sht_u(:,1,iii),'b-');
hold on
plot3(USht_x(:,1,iii),USht_x(:,2,iii),USht_u(:,1,iii),'k-');
xlabel('X');
ylabel('Y');
zlabel('U');
grid on;


N = size(Sht_x, 3);
figure(4);
for iii = 1:N
    xvec = repmat(SOM1(iii), size(Sht_x,1), 1);
    hold on;
    plot3(xvec, Sht_x(:,1,iii), Sht_x(:,2,iii), 'b-');
    plot3(xvec, USht_x(:,1,iii), USht_x(:,2,iii), 'k-');
end
xlabel('k_{OM1}');
ylabel('X');
zlabel('Y');
grid on;
view(3);

St_x   = coco_bd_col(bd22, 'SX.T_x');
USt_x  = coco_bd_col(bd22, 'USX.T_x');
St_u   = coco_bd_col(bd22, 'SX.T_u');
USt_u  = coco_bd_col(bd22, 'USX.T_u');

iii = 1;
figure(5)
plot3(St_x(:,1,iii), St_x(:,2,iii), St_u(:,1,iii), 'b-');
hold on
plot3(USt_x(:,1,iii), USt_x(:,2,iii), USt_u(:,1,iii), 'k-');
xlabel('X');
ylabel('Y');
zlabel('U');
grid on;

N = size(St_x, 3);
figure(6);
for iii = 1:N
    xvec = repmat(SOM1(iii), size(St_x,1), 1);
    plot3(xvec, St_x(:,1,iii), St_x(:,2,iii), 'b-');
    hold on;
    plot3(xvec, USt_x(:,1,iii), USt_x(:,2,iii), 'k-');
end
xlabel('k_{OM1}');
ylabel('X');
zlabel('Y');
grid on;
view(3);

function [data y J] = my_DC(prob,data,u)
SC  = u(1:end/2);
USC = u(end/2+1:end);
y = abs(norm(USC) - norm(SC));
J = [];
end
