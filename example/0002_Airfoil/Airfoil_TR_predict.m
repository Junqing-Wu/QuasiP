function [Xk_MF,OM,y_qp] = Airfoil_TR_predict(Psi,Om1,P1,eps,y_po)
% AIRFOIL_TR_PREDICT  Standalone copy of the TR-to-2-torus predictor.
% Copied from PO_TR2QP_d2 and its Fourier helpers in ode_TR2QuasiP.
% This example-local function does not call or modify ode_TR2QuasiP.
% Reference: QuasiP article, Section 5.2 and Eq. (35).
%
% Inputs: Psi - state-transition history along the periodic orbit;
%         Om1 - periodic base frequency; P1 - periodic multiplication factor;
%         eps - perturbation amplitude; y_po - rows of [q, qdot].
% Outputs: Xk_MF - predicted Fourier representation; OM - two base frequencies;
%          y_qp - full samples of the predicted torus (an initial guess only).
% The prediction, multiplier selection, and Fourier conversion are retained
% from the toolbox implementation so the original visualization is reproducible.
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
    Xk_MF = QP2X_MF(y_qp,Y,P1);
end
end

function Xk_MF = QP2X_MF(y_qp,Y,P1)
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

[XC_ii,A_ii] = sort(X2_C(:,1),'descend');
idx  = abs(XC_ii) > max(abs(XC_ii))*1e-7;
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