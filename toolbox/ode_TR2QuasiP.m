function prob = ode_TR2QuasiP(prob, oid, varargin)
% ODE_TR2QUASIP   Start a QuasiP branch from a torus bifurcation.
%
% PROB = ODE_TR2QUASIP(PROB, OID, RUN, [SOID], LAB, [EPS], [DFDTH])
%
% Read a torus bifurcation of a periodic-orbit branch, perturb its critical
% eigenvector to create an initial two-frequency torus, and construct a new
% QuasiP continuation problem from that initial guess.
%
% On input:
%
% PROB  : COCO continuation problem structure.
% OID   : Target QuasiP object identifier.
% RUN   : Source periodic-orbit continuation run.
% SOID  : Optional source object identifier.
% LAB   : Label of the torus bifurcation point.
% EPS   : Initial perturbation amplitude; default is 1e-3.
% DFDTH : Optional excitation-phase derivative callback.
%
% On output:
%
% PROB  : COCO problem containing the new QuasiP branch.
%
% The torus-point data are converted to Xk_MF and Om before being passed to
% ODE_ISOL2QUASIP.
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
[sol,data] = QuasiP_read_solution(args.soid,args.run,args.lab);
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
bd     = get_TRbd(args);
Psi    = coco_bd_col(bd, 'Psi');
Om1    = coco_bd_col(bd, 'Om');
n      = data.ode_opts.n;
x_po   = coco_bd_col(bd, 'HT_x');
fullOutput = ~isempty(x_po) && size(x_po,2)==n && isequal(data.amp_opts.outdof(:),(1:n)');
if data.ODE == 2 && fullOutput && data.amp_opts.isuamp
    u_po = coco_bd_col(bd, 'HT_u');
    fullOutput = ~isempty(u_po) && isequal(size(u_po),size(x_po));
end
if data.ODE==2 && ~data.amp_opts.isuamp
    fullOutput = false;
end
if ~fullOutput && strcmp(data.PHI.param,'MP')
    % Saved plots may contain only selected DOFs or no velocity history.
    Amp = PostP_amp(data,sol.u(data.C_idx),Om1,1:n,data.ODE==2);
    x_po = Amp.HT_x;
    if data.ODE == 2
        u_po = Amp.HT_u;
    end
elseif ~fullOutput || (data.ODE==2 && ~data.amp_opts.isuamp)
    error('ode_TR2QuasiP:IncompleteOrbit','Save full periodic states and velocities before switching from this TR.');
end
if data.ODE == 2
    y_po = [x_po,u_po];
else
    y_po = x_po;
end
[Xk_MF0,OM0] = PO_TR2QP_d2(Psi,Om1,data.PHI.P{1},eps,y_po,data.ode_opts.e);
if data.ODE == 2
   Xk_MF0.Xc = Xk_MF0.Xc(:,1:end-data.ode_opts.n);
end
k12 = Xk_MF0.k_wave(:,1:2);
mask = ismember(k12, [1 -1; 1 1], 'rows');
if data.ode_opts.e == 0 && ~any(mask)
    Xk_MF0.k_wave(:,1:2) = Xk_MF0.k_wave(:,1:2)/k12(1:2,:);
    OM0          = k12(1:2,:)*OM0;
end
%--------------------------------------------------------------------------
pname     = data.pnames_np;
p0        = coco_bd_col(bd, 'p');
%--------------------------------------------------------------------------
prob      = ode_isol2QuasiP(prob, oid, funcs{:}, Xk_MF0, OM0, pname, p0);
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

function [Xk_MF,OM,y_qp]  = PO_TR2QP_d2(Psi,Om1,P1,eps,y_po,e)
[N,L] = size(Psi,1,3);
S     = size(y_po,1);
m     = S/(L-1);
mm    = (L-1)/S;
assert(mod(m,1)==0 || mod(mm,1)==0, ...
    'the ratio of S_wave in PHI and L in sta_opts is not an integer.');
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
    Xk_MF = QP2X_MF(y_qp,Y,P1,e);
end
end

function Xk_MF = QP2X_MF(y_qp,Y,P1,e)
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
k_wave = kQ';

if e > 0
    % A forced beam's first displacement can have no critical sidebands.
    [XC_ii,A_ii] = sort(max(abs(X2_C),[],2),'descend');
else
    [XC_ii,A_ii] = sort(X2_C(:,1),'descend');
end
idx = abs(XC_ii) > max(abs(XC_ii))*1e-7;
A = A_ii(idx);

Xk_MF.k_wave = k_wave(A,:);
Xk_MF.Xc     = X2_C(A,:);
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
