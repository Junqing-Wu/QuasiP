% DEMO_TR2QP  Periodic and quasi-periodic rotor-stator rubbing branches.
% Model: qddot + 0.1*qdot + 0.04*q = F_contact(q,qdot,mu,OM)
%        + OM^2*[cos(theta1);sin(theta1)], q=[x;y], theta1=omega1*t.
% The glue constraint enforces omega1=OM, with p=[mu;OM].
% Start from integration at omega1=0.2, continue periodic orbits as qp1,
% then initialize a two-frequency torus at the selected TR point as qp2.
% Save k_OM1=0.30,0.35,0.40,0.45 events on both branches.
% After both runs finish, execute figure_RS_events for phase portraits and spectra.
% Event availability depends on actual continuation coverage; no event is fabricated.
% The stability settings below are retained from the original example and are costly.
% Initialize COCO/QuasiP first and work in this example directory.
% -------------------------------------------------------------------------
clc;
clear;
% Physical parameters and prescribed excitation frequency.
rho = 1;
Om1 = 0.2;
OM0 = rho*Om1;
pname = {'mu','OM'};
p0 = [0.2;OM0(1)];
%%
% initial condition
% Integrate the physical model and retain its late-time spectrum for initialization.
x0 = zeros(4,1)+0.0001;
Num_period = 2^9;
T1 = 2*pi/(OM0(1));
S_period = 2^10;
dt = T1/S_period;
bbb = 2;
tspan = [0 Num_period*T1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
isfwithu1 = true; [t1,x1] = ode15s(@(t,x) my_dxdt(t,x,p0,OM0,@RS1_F,@RS1_L,isfwithu1), tspan, x0, opts);
tt = t1; xx = x1;
t = Num_period/bbb*T1:dt:Num_period*T1-dt;
x = zeros(length(t), size(x0,1)); 
for j = 1:size(x0,1)
    x(:,j) = interp1(tt, xx(:,j), t, 'spline'); 
end
X = fft(x)/length(x(:,1)); X(2:end,:)=X(2:end,:)*2;
na = length(x(:,1)); ff = (0:na-1)'/na*2*pi/dt;
[pks,locs] = findpeaks(log(abs(X(1:end/2,1))));
[pks_sorted, idx_sorted] = sort(pks, 'descend');
keep_idx = pks_sorted >= (pks_sorted(1) - 8);  
pks_final = pks_sorted(keep_idx);
loc_final = locs(idx_sorted(keep_idx));
ff_final = ff(loc_final,:)/(OM0(1));
X_final  = X(loc_final,:);
rho = ff_final(1);
Xk_MF = solve_HOrder(X_final,ff_final,rho,5,1);
Xk_MF2 = Xk_MF;
Xk_MF2.Xc = Xk_MF.Xc(:,1:(3-size(Xk_MF.k_wave,2)));
% Convert the Fourier initial guess into a periodic m-VCF continuation problem.

OM0 = rho*Om1;
prob1 = coco_prob();
ode_opts = struct('n',2, 'd',1, 'e',1, 'isfwithu',true);
prob1 = coco_set(prob1,'QuasiP','ode_opts',ode_opts); 
prob1 = coco_set(prob1,'cont','PtMX',3000);
prob1 = coco_set(prob1,'cont','h_min',0.01,'h_max',0.05);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
sta_opts = struct('issta',true,'S',1,'L1',2^9,'N_Ly',[],'Omt',1);
prob1 = coco_set(prob1,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 

PHI_INITIAL = struct( 'param', 'MP', ...
    'method',  'mVCF'           , ...     
    'type',    { {'HB'} }, ...   
    'k',       { {5} }, ...        
    'P',       { {1} }, ...    
    'S',       { {2^8} }, ...
    'S_wave',  { {2^8} }, ...
    'OmOd',    1);
prob1 = coco_set(prob1, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
funcs = {@RS1_F, @RS1_DFDX, @RS1_DFDP, @RS1_L};

prob1 = ode_isol2QuasiP(prob1, '', funcs{:}, Xk_MF2, OM0, pname , p0);
[data, uidx] = coco_get_func_data(prob1, 'QuasiP', 'data', 'uidx');
prob1 = coco_add_glue(prob1,'OMp2',data.kOM_idx,data.p_idx(2));
% Localize and save the requested frequency sections on this branch.
% These event labels are consumed by figure_RS_events.m.
event_values = [0.30,0.35,0.40,0.45];
event_names = {'Om1=0.30','Om1=0.35','Om1=0.40','Om1=0.45'};
for ievent = 1:numel(event_values)
    prob1 = coco_add_event(prob1,event_names{ievent},'k_OM1',event_values(ievent));
end
bd1 = coco(prob1, 'qp1', [], 1, { 'k_OM1' 'OM' 'mu' 'QuasiP.test.TR'}, [0.15 0.5]);
% close(11111)

%% Initialize the quasi-periodic branch from a saved periodic TR point
lab  = coco_bd_labs(bd1, 'TR');
assert(~isempty(lab), 'No saved TR point in qp1; inspect the periodic continuation log.');
lab  = max(lab);
prob2 = coco_prob();
prob2 = coco_set(prob2,'corr','ItMX',50);
amp_opts    = struct('isamp',true,'isuamp',true,'isht',true,'outdof',[1,2],'Nperiod',1,'Nsampling',2^8);
sta_opts = struct('issta',true,'S',2^10,'L1',2^10,'N_Ly',10000,'Omt',1);
prob2 = coco_set(prob2,'QuasiP','amp_opts',amp_opts,'sta_opts',sta_opts); 
PHI_INITIAL = struct( 'param', 'MP', ...
    'method',   'mVCF'           , ...     
    'type',    { {'HB','HB'} }, ...  
    'k',       { {11,11} }      , ...        
    'P',       { {1,1} }      , ...    
    'S',       { {2^8,2^8} }  , ...
    'S_wave',  { {2^8,2^8} }  , ...
    'OmOd' ,   [1,2]);
prob2 = coco_set(prob2, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 
prob2 = ode_TR2QuasiP(prob2, '', 'qp1', '', lab);
[data, uidx] = coco_get_func_data(prob2, 'QuasiP', 'data', 'uidx');
prob2 = coco_add_glue(prob2,'OMp2',data.kOM_idx,data.p_idx(2));
% Localize and save the requested frequency sections on this branch.
% These event labels are consumed by figure_RS_events.m.
event_values = [0.30,0.35,0.40,0.45];
event_names = {'Om1=0.30','Om1=0.35','Om1=0.40','Om1=0.45'};
for ievent = 1:numel(event_values)
    prob2 = coco_add_event(prob2,event_names{ievent},'k_OM1',event_values(ievent));
end
coco(prob2, 'qp2', [], 1, { 'k_OM1' 'u_om2' 'OM' 'mu'}, [0.28927 1]);
% close(11111)


