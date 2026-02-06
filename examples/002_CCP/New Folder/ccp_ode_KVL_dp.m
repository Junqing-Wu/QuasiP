function df_dp = ccp_ode_KVL_dp(t, x, p)

    Q1  = x(1,:);   % Q_loop1
    Q2  = x(2,:);   % Q_loop2
    Q3  = x(3,:);   % Q_loop3
    I2  = x(4,:);   % I_loop2
    I3  = x(5,:);   % I_loop3
    Qs1 = x(6,:);   % Q_s1
    Qs2 = x(7,:);   % Q_s2

    Cm1 = p(1,:);
    Cm2 = p(2,:);
    wrf = p(3,:);

    par.Rrf   = 50;
    par.Lm2   = 1500e-9;
    par.Vrf   = 100;         % [V]   RF source voltage amplitude (25 W)
    
    % ---------------- 解包状态变量 ----------------
    Q1  = x(1,:);   % Q_loop1
    Q2  = x(2,:);   % Q_loop2
    
    % 初始化导数向量
    N       = numel(Cm1);
    df_dp   = zeros(7,2,N);
    
    df_dp(1,1,:) = (Q1 - Q2) ./ (Cm1.^2 .* par.Rrf);
    df_dp(4,1,:) = (Q2 - Q1) ./ (Cm1.^2 .* par.Lm2);
    
    df_dp(4,2,:) = Q2./(Cm2.^2 .* par.Lm2);

    df_dp(1,3,:) = (par.Vrf .* t .* sin(wrf .* t)) ./ par.Rrf;
    % 其余分量（dQ2/dQ3/dI3/dQs1/dQs2）与Cm1无关，导数为0
end