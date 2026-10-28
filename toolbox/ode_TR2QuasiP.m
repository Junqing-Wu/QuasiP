function prob = ode_TR2QuasiP(prob, oid, varargin)
%ODE_TR2QuasiP   Switch to branch of QuasiP at torus (TR) bifurcation.
%
% PROB = ODE_TR2QUASIP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [EPS [DFDTH]] }

% Start a continuation of quasiperiodic orbits emanating from a torus bifurcation
% that was obtained and saved to dis in a previous continuation along a
% branch of periodic orbits with name RUN. To start from a saved TR
% bifurcation, at least the name RUN of the continuation run and the
% solution label LAB must be given. The label LAB must be the label of a
% TR bifurcation points.

% PROB   - Continuation problem structure.
% OID    - Target object instance identifier (string).
% RUN    - Run identifier (string or cell-array of strings).
% SOID   - Source object instance identifier (string, optional).
% LAB    - Solution label (integer).
% EPS    - Amount of peturbation, 1e-3 by default
% @DFDTH - Optional function handle to Jacobian w.r.t. theta

grammar   = 'RUN [SOID] LAB [EPS] [DFDTH]';
args_spec = {
     'RUN', 'cell', '{str}',       'run',   {},   'read',  {}
    'SOID',     '',   'str',      'soid',  oid,   'read',  {}
     'LAB',     '',   'num',       'lab',   [],   'read',  {}
     'EPS',     '',   'num',       'eps',   [],   'read',  {}
     'DFDTH',   '',  '@|[]',  'dfdthhan',   [],   'read',  {} 
     };
opts_spec = [];
[args, ~] = coco_parse(grammar, args_spec, opts_spec, varargin{:});
if isempty(args.eps)
    eps = 1e-3;
else
    eps = args.eps;
end
fth = args.dfdthhan;
%--------------------------------------------------------------------------
prob0  = coco_prob();
prob0  = ode_QuasiP2QuasiP(prob0, '' , args.run, args.soid, args.lab);
data   = coco_get_func_data(prob0, 'QuasiP', 'data'); clear prob0;
%--------------------------------------------------------------------------
ode_opts  = data.ode_opts; ode_opts.d = 2;
prob      = coco_set(prob,'QuasiP','ode_opts',ode_opts);
%--------------------------------------------------------------------------
if isempty(fth)
    funcs = {data.f, data.fx, data.fp, data.l, data.fth};
else
    funcs = {data.f, data.fx, data.fp, data.l, fth};
end
%--------------------------------------------------------------------------
if data.ODE == 2
   assert(data.amp_opts.isuamp, ...
  'run=%s data.amp_opts.isuamp must be true', args.run{1});
end
bd     = get_TRbd(args);
Psi    = coco_bd_col(bd, 'Psi');
Om1    = coco_bd_col(bd, 'Om');
x_po   = coco_bd_col(bd, 'HT_x');
if data.ODE == 2
    u_po   = coco_bd_col(bd, 'HT_u');
    y_po = [x_po,u_po];
else
    y_po = x_po;
end
[X_MF0,OM0]  = PO_TR2QP_d2(Psi,Om1,data.PHI.P{1},eps,y_po);

% %-----------------------
% figure(1)
% plot3([y_po(:,1);y_po(1,1)],[y_po(:,4);y_po(1,4)],[y_po(:,2);y_po(1,2)],'k--','LineWidth',1.5)
% hold on
% [X_MF0,OM0,y_qp]  = PO_TR2QP_d2(Psi,Om1,data.PHI.P{1},0.05,y_po);
% y_qp = reshape(y_qp,256,256,[]);
% y_qp = y_qp(1:8:end,1:8:end,:);
% y_qp = reshape(y_qp,32*32,[]);
% plot3(y_qp(:,1),y_qp(:,4),y_qp(:,2),'g.')
% 
% [X_MF0,OM0,y_qp]  = PO_TR2QP_d2(Psi,Om1,data.PHI.P{1},0.1,y_po);
% y_qp = reshape(y_qp,256,256,[]);
% y_qp = y_qp(1:8:end,1:8:end,:);
% y_qp = reshape(y_qp,32*32,[]);
% plot3(y_qp(:,1),y_qp(:,4),y_qp(:,2),'r.')
% 
% [X_MF0,OM0,y_qp]  = PO_TR2QP_d2(Psi,Om1,data.PHI.P{1},0.2,y_po);
% y_qp = reshape(y_qp,256,256,[]);
% y_qp = y_qp(1:8:end,1:8:end,:);
% y_qp = reshape(y_qp,32*32,[]);
% plot3(y_qp(:,1),y_qp(:,4),y_qp(:,2),'b.')
% 
% xlabel('$q_1$','Interpreter','latex');
% ylabel('$\dot{q}_1$','Interpreter','latex');
% zlabel('$q_2$','Interpreter','latex');
% legend('PO-unstable','QP-$\varepsilon=0.05$','QP-$\varepsilon=0.10$','QP-$\varepsilon=0.20$','Location','best');
% grid on;
% plot_setting_paper()
% set(findall(gcf,'type','line'),'linewidth',3)
% legend boxon;
% %---------------------------------------

if data.ODE == 2
   X_MF0 = X_MF0(:,1:end-data.ode_opts.n);
end
k12 = X_MF0(:,1:2);
mask = ismember(k12, [1 -1; 1 1], 'rows');
if isempty(find(mask, 1))
    X_MF0(:,1:2) = X_MF0(:,1:2)/k12(1:2,:);
    OM0          = k12(1:2,:)*OM0;
end
%--------------------------------------------------------------------------
pname     = data.pnames_np;
p0        = coco_bd_col(bd, 'p');
%--------------------------------------------------------------------------
prob      = ode_isol2QuasiP(prob, oid, funcs{:}, X_MF0, OM0, pname, p0);
end


function TRbd = get_TRbd(args)
bd   = coco_bd_read(args.run);
lab  = coco_bd_col(bd, 'LAB');
mask = cellfun(@(x) ~isempty(x) && (x == args.lab), lab);
pos  = find(mask)+1;
assert(~isempty(pos), ...
  'get_TRbd: run=%s no TR LAB=%d', args.run{1}, args.lab);
TRbd = bd([1,pos],:);
type  = coco_bd_col(TRbd, 'TYPE');
assert(strcmp(type{1},'TR'), ...
  'get_TRbd: run=%s, LAB=%d, TYPE = %s, must be TR', args.run{1}, args.lab, type{1});
end

function [X_MF,OM,y_qp]  = PO_TR2QP_d2(Psi,Om1,P1,eps,y_po)
[N,L] = size(Psi,1,3);
S     = size(y_po,1);
m     = S/(L-1);
mm    = (L-1)/S;
assert(mod(m,1)==0 || mod(mm,1)==0, ...
    'the ratio of S_wave in PHI and L in stable_opts is not an integer.');
M = Psi(:,:,end);
[vector,value] = eig(M);
value  = diag(value);
idx    = abs(real(value)-1) < 1e-4; value  = value(~idx); vector = vector(:,~idx);
idx    = abs(real(value)+1) < 1e-4; value  = value(~idx); vector = vector(:,~idx);
idx    = abs(abs(value)-1)  < 1e-4; psi = value(idx); vec = vector(:,idx);
alpha  = atan(imag(psi)./real(psi));
idx    = alpha > 0; psi = psi(idx); v0 = vec(:,idx); alpha = alpha(idx);
if m == 1 
    Psi  = Psi(:,:,1:end-1);  Y = S/P1;
elseif m > 1
   y_po = y_po(1:m:end,:); Psi  = Psi(:,:,1:end-1); Y = (L-1)/P1;
elseif m < 1
   Psi  = Psi(:,:,1:mm:end-1);Y = S/P1;
end
if numel(psi) ~= 1
    warning('Not a pair of conjugate eigenvalues on the unit circle in TR2QuasiP.');
else
    OM  =  [Om1;alpha*Om1/(2*P1*pi)];
    v   = zeros(N,Y*P1);
    for ii = 1:Y*P1
        v(:,ii) = exp(-1i*alpha*(ii-1)/(Y*P1))*Psi(:,:,ii)*v0;
    end
    rv = real(v);
    iv = imag(v);
    y1_bar   = kron(y_po',[1;0;0]) + kron(eps * rv,[0;1;0])  + kron(- eps * iv,[0;0;1]) ;
    Tau2     = 0:2*pi/(Y):2*pi-2*pi/(Y);
    Phi2_bar = [ones(Y,1),cos(Tau2)',sin(Tau2)'];
    y_qp = reshape( transpose( kron(eye(N),Phi2_bar)*y1_bar ), [], N ) ;
    X_MF = QP2X_MF(y_qp,Y,P1);
end
end

function X_MF = QP2X_MF(y_qp,Y,P1)
k1 = 1:1/P1:7;
k = {k1,1};
N = {numel(k1),1};
kQ = Com_kQ(k,N);
U = size(kQ,2)*2+1;

tau_1im1 = (0:2*pi/Y/P1:2*pi-2*pi/Y/P1)';
S_1im1 = Y*P1;
for ii = 2:2
    tau_i = (0:2*pi/Y:2*pi-2*pi/Y)';
    tau_1im1 = [ kron( ones(Y,1),tau_1im1 ) , kron( tau_i,ones(S_1im1,1) ) ];
    S_1im1 = S_1im1*Y;
end
tau = tau_1im1;

GAM      = zeros( Y^2*P1, U );
GAM(:,1) = ones( Y^2*P1,1 );
GAM(:,2:2:end) = cos(tau*kQ);
GAM(:,3:2:end) = sin(tau*kQ);
iGAM = GAM';
iGAM(2:end,:) = 2*iGAM(2:end,:);
iGAM = iGAM/(Y^2*P1);

X2 = iGAM*y_qp;
X2_C = [X2(1,:) ; X2(2:2:end,:) - 1i*X2(3:2:end,:)];
kQ = [[0;0],kQ];
X_MF = [kQ',X2_C];

[XC_ii,A_ii] = sort(X2_C(:,1),'descend');
idx  = abs(XC_ii) > max(abs(XC_ii))*1e-7;
A = A_ii(idx);

X_MF = X_MF(A,:);
end

function kQ = Com_kQ(k,N)
kQ = [0,k{1}];
d = length(N);
for ii = 2:d
    kQ = kQ(:,2:end);
    NQ = size(kQ,2);
    OO = zeros(ii-1,1);
    kQ  = [ [kron(ones(1,N{ii}+1),[OO,kQ]);kron([0,k{ii}],ones(1,NQ+1))],...
          [kron(ones(1,N{ii}),kQ)  ;kron(-k{ii},ones(1,NQ))] ];
end
kQ = kQ(:,2:end);
end
