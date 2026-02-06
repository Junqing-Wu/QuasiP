function dxdt = ccp_ode_dual(t, x, params)
    % CCP_ODE_DUAL 双频动力学方程组的ODE函数（向量化计算）
    % 输入：
    %   t       - 时间标量/向量（s），ODE求解器自动传入
    %   x       - 状态列向量（9×1），对应状态量：
    %             x(1)=Q_loop1, x(2)=Q_loop2, x(3)=Q_loop3, x(4)=Q_loop4,
    %             x(5)=I_loop2, x(6)=I_loop3, x(7)=I_loop4,
    %             x(8)=Q_s1,    x(9)=Q_s2
    %   params  - 参数字符体，包含所有方程参数（见下方参数列表）
    % 输出：
    %   dydt    - 状态导数列向量（9×1），与y一一对应
    %
    % 参数字符体params需包含的字段：
    %   VsHF_star, omega_HF,  % 高频激励：幅值、角频率
    %   VsLF_star, omega_LF,  % 低频激励：幅值、角频率
    %   Cm1HF, RsHF, RmHF, Cm2HF, LmHF,  % 高频支路参数
    %   Cm1LF, Cm2LF, RmLF, RsLF, LmLF,  % 低频主支路参数
    %   Rpl, Rf, Lpl, Lf,     % 等离子体/滤波支路参数
    %   Cs10, Cs20, Te,       % 鞘层电容/电子温度参数
    %   Ii1, Ie10,            % 鞘层1离子/电子电流参数
    %   Ii2, Ie20             % 鞘层2离子/电子电流参数
    
    %% 1. 解包状态向量（提高代码可读性，避免下标混乱）
    Q_loop1 = x(1,:); 
    Q_loop2 = x(2,:); 
    Q_loop3 = x(3,:); 
    Q_loop4 = x(4,:);
    I_loop2 = x(5,:); 
    I_loop3 = x(6,:); 
    I_loop4 = x(7,:);
    Q_s1    = x(8,:); 
    Q_s2    = x(9,:);
    
    %% 2. 解包参数字符体
    % 高频激励参数
    VsHF_star = params.VsHF_star; omega_HF = params.omega_HF;
    % 低频激励参数
    VsLF_star = params.VsLF_star; omega_LF = params.omega_LF;
    % 高频支路参数
    Cm1HF = params.Cm1HF; RsHF = params.RsHF; RmHF = params.RmHF;
    Cm2HF = params.Cm2HF; LmHF = params.LmHF;
    % 低频主支路参数
    Cm1LF = params.Cm1LF; Cm2LF = params.Cm2LF; RmLF = params.RmLF;
    RsLF  = params.RsLF;  LmLF  = params.LmLF;
    % 等离子体/滤波支路参数
    Rpl = params.Rpl; Rf = params.Rf; Lpl = params.Lpl; Lf = params.Lf;
    % 鞘层与电子参数
    Cs10 = params.Cs10; Cs20 = params.Cs20; Te = params.Te;
    % 鞘层电流参数
    Ii1 = params.Ii1; Ie10 = params.Ie10;
    Ii2 = params.Ii2; Ie20 = params.Ie20;
    
    %% 3. 计算时变激励电压（向量化：点乘兼容t为向量）
    VsHF = VsHF_star .* sin(omega_HF .* t);
    VsLF = VsLF_star .* sin(omega_LF .* t);
    
    %% 5. 逐行计算状态导数（严格遵循原方程，括号匹配，向量化点运算）
    % 5.1 dQ_loop1/dt
    dQ_loop1_dt = ( (Q_loop4 - Q_loop1)./Cm1HF - VsHF ) ./ RsHF;
    
    dQ_loop2_dt = I_loop2;  % dQ_loop2/dt = I_loop2
    
    dQ_loop3_dt = I_loop3;  % dQ_loop3/dt = I_loop3
    
    dQ_loop4_dt = I_loop4;  % dQ_loop4/dt = I_loop4
    
    % 5.5 dI_loop2/dt
    term_I2 = -(RmLF + RsLF).*dQ_loop2_dt ...
              - (1./Cm1LF + 1./Cm2LF).*Q_loop2 ...
              + (Q_loop3 - Q_loop4)./Cm1LF ...
              - VsLF;
    dI_loop2_dt = term_I2 ./ LmLF;
    
    % 5.6 dI_loop3/dt
    term_I3 = -(Rpl + Rf).*dQ_loop3_dt ...
              + (Q_loop2 - Q_loop3 + Q_loop4)./Cm1LF ...
              - Q_s1.^2 ./ (2.*Cs10.^2 .* Te) ...
              + Q_s2.^2 ./ (2.*Cs20.^2 .* Te);
    dI_loop3_dt = term_I3 ./ (Lpl + Lf);
    
    % 5.7 dI_loop4/dt
    term_I4 = -RmHF.*dQ_loop4_dt ...
              + Q_loop1./Cm1HF ...
              + (Q_loop3 - Q_loop2)./Cm1LF ...
              - (1./Cm1LF + 1./Cm2HF + 1./Cm1HF) .* Q_loop4;
    dI_loop4_dt = term_I4 ./ LmHF;
    
    % 5.8 dQ_s1/dt（指数项：鞘层电子电流的玻尔兹曼分布）
    exp_s1 = exp( - (Q_s1 ./ (sqrt(2).*Cs10.*Te)) .^2 );
    dQ_s1_dt = dQ_loop3_dt - Ii1 + Ie10 .* exp_s1;
    
    % 5.9 dQ_s2/dt（注意dQ_loop3/dt前的负号）
    exp_s2 = exp( - (Q_s2 ./ (sqrt(2).*Cs20.*Te)) .^2 );
    dQ_s2_dt = -dQ_loop3_dt - Ii2 + Ie20 .* exp_s2;
    
    %% 6. 转换为列向量（MATLAB ODE求解器要求导数为列向量）
    dxdt = [dQ_loop1_dt;...
            dQ_loop2_dt;...
            dQ_loop3_dt;...
            dQ_loop4_dt;...
            dI_loop2_dt;...
            dI_loop3_dt;...
            dI_loop4_dt;...
            dQ_s1_dt;...
            dQ_s2_dt];

end