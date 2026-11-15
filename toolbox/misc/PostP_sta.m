function Sta = PostP_sta(data,Omii,pii,Amp)
%POSTP_STA  Dispatch stability processing and classify the current solution.
%
% Select the MP or MT transition-data constructor, then compute characteristic
% multipliers for d=1 or Lyapunov quantities for d>1. The result is used by
% QuasiP_d1_test and by the COCO branch-data slot.
%
% On input:
%   data    : QuasiP model, PHI, and stability settings.
%   Omii    : Current base-frequency vector.
%   pii     : Current model-parameter vector.
%   Amp     : Amplitude/reconstruction data.
%
% On output:
%   Sta     : Stability data, transition matrices, multipliers/exponents, and labels.

switch data.PHI.param
    case{'MP'}
        Psi = MP_PostP_sta_Psi(data,Omii,pii,Amp);
    case{'MT'}
        switch data.PHI.method
            case{'DM'}
            case{'SHT'}
                Psi = MT_SHT_PostP_sta_Psi(data,Omii,pii,Amp);
        end
end

Perii   = cell2mat(data.PHI.P);
if data.ode_opts.d == 1
    Sta.Psi = Psi;
    [~,eigvalue] = eigs(Psi(:,:,end),size(Psi,2));
    eigvalue = diag(eigvalue);

%     T = 2*pi /Omii;
%     R = 1/T*logm(Psi(:,:,end))
%     [~,eigvalueA] = eigs(R,size(Psi,2));
%     eigvalueA = diag(eigvalueA);

    Flo_mul = eigvalue;
    if max( abs(Flo_mul) )> 1
        Sta.Stab = 0;
    else
        Sta.Stab = 1;
    end
    Flo_exp = log(Flo_mul)*Omii/(Perii*2*pi);

%     [Flo_exp,eigvalueA]

    Sta.Ly = [Flo_mul,Flo_exp];
else
    Sta.Psi = Psi(2,:);
    Omi      = data.sta_opts.Omt;
    NumOmi   = numel(Omi);
    Sta.Ly   = [];
    Sta.Stab = [];
    for iii = 1:NumOmi
        Om_t    = Omi(iii);
        Psi_iii = Psi{1,iii};
        N_Ly     = data.sta_opts.N_Ly;
        N        = sqrt(size(Psi_iii,2));
        E        = eye(N);
        lnR      = zeros(1,N);
        Lyapunov = zeros(N_Ly,N);
        dTi = 2*Perii(Om_t)*pi/Omii(Om_t);
        for ii = 1:N_Ly
            Mii = reshape(Psi_iii(ii,:),N,N);
            %         Mii = pinv(Mii);
            %         Mii = Mii';
            E = Mii*E;
            [E, R] = householder_qr(E);

            lnR = lnR + log(abs(diag(R)'));
            if any(isnan(R), 'all') || any(isinf(R), 'all')
                error('QR factorization returned NaN or Inf.');
            end
            Lyapunov(ii,:) = lnR / ii / dTi;
        end
        Ly     = sum(Lyapunov(end-20:end,:))/21;
        Sta.Ly = [Sta.Ly,Ly];
        if max( Ly ) > 0.001
            Sta.Stab = [Sta.Stab,0];
        else
            Sta.Stab = [Sta.Stab,1];
        end
    end
end
end


function [Q, R] = householder_qr(A)
% Input:
% A: Matrix to factorize.
% Output:
% Q: Orthogonal factor.
% R: Upper-triangular factor.

[m, n] = size(A);   % Get the matrix dimensions.
Q = eye(m);         % Initialize Q as the identity.
R = A;              % Initialize R with A.

% Apply Householder reflections column by column.
for k = 1:min(m, n)
    % Extract the active part of column k.
    x = R(k:m, k);  % Extract the active part of column k.

    % Construct the Householder reflection vector v.
    e1 = zeros(length(x), 1);
    e1(1) = 1;  % Construct the first basis vector e1.
    v = sign(x(1)) * norm(x) * e1 + x;  % Construct the reflection vector.
    v = v / norm(v);  % Normalize the reflection vector.

    % Construct the Householder matrix H_k.
    H_k = eye(m - k + 1) - 2 * (v * v');  % Form the Householder matrix.

    % Apply H_k to R.
    R(k:m, :) = H_k * R(k:m, :);  % Update the active rows of R.

    % Construct the full Householder matrix H_full.
    H_full = eye(m);
    H_full(k:m, k:m) = H_k;

    % Update Q.
    Q = Q * H_full;
end

end

