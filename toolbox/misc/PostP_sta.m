function Sta = PostP_sta(data,Omii,pii,Amp)

switch data.PHI.parametrization_tori
    case{'MP'}
        Psi = MP_PostP_sta_Psi(data,Omii,pii,Amp);
    case{'MT'}
        switch data.PHI.type
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
    Omi      = data.stable_opts.Omi;
    NumOmi   = numel(Omi);
    Sta.Ly   = [];
    Sta.Stab = [];
    for iii = 1:NumOmi
        Om_t    = Omi(iii);
        Psi_iii = Psi{1,iii};
        N_Ly     = data.stable_opts.N_Ly;
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
                error('QR分解结果包含NaN或无穷值！');
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
% 输入:
% A: 需要进行 QR 分解的矩阵
%
% 输出:
% Q: 正交矩阵
% R: 上三角矩阵

[m, n] = size(A);   % 获取矩阵 A 的大小
Q = eye(m);         % 初始化 Q 为单位矩阵
R = A;              % 初始化 R 为矩阵 A

% 逐列进行 Householder 变换
for k = 1:min(m, n)
    % 提取 A 的第 k 列向量的下部分
    x = R(k:m, k);  % 提取当前列中从第 k 行开始的部分

    % 构造 Householder 反射向量 v
    e1 = zeros(length(x), 1);
    e1(1) = 1;  % 构造标准基向量 e1
    v = sign(x(1)) * norm(x) * e1 + x;  % 构造反射向量
    v = v / norm(v);  % 归一化

    % 构造 Householder 矩阵 H_k
    H_k = eye(m - k + 1) - 2 * (v * v');  % 计算 Householder 变换矩阵

    % 将 H_k 应用于矩阵 R
    R(k:m, :) = H_k * R(k:m, :);  % 更新 R 中从 k 行开始的部分

    % 构造完整的 Householder 变换矩阵 H_full
    H_full = eye(m);
    H_full(k:m, k:m) = H_k;

    % 更新 Q 矩阵
    Q = Q * H_full;
end

end


