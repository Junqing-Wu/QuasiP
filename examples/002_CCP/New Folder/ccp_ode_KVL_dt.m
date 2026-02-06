function dt_dx = ccp_ode_KVL_dt(t, x, p)

    par.Rrf   = 50;          % [Ohm] Equivalent source internal resistance
    par.Vrf   = 100;         % [V]   RF source voltage amplitude (25 W)
    par.frf   = 13.56e6;     % [Hz]  RF source frequency
    % par.wrf   = 2*pi*par.frf;% [rad/s]
    
    wrf = p(3,:);
    % 初始化导数向量
    N     = numel(t);
    dt_dx = zeros(7, N);
    
    % 计算vsrc对t的导数（仅dQ1含时间项）
    % dQ1对t的偏导：-d(vsrc)/dt / Rrf
    dt_dx(1,:) = (par.Vrf .* wrf .* sin(wrf .* t)) / par.Rrf;
    
    % 其余分量（dQ2/dQ3/dI2/dI3/dQs1/dQs2）与t无关，导数为0
end