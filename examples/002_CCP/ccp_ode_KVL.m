function dx = ccp_ode_KVL(t, x, p)
    %CCP_ODE  dx/dt = f(t,x) for CCP equivalent circuit (charge-state form)
    %% Params
    par.Cs10 = 1.08e-10;     % [F]  RF sheath S1 equivalent capacitance constant term
    par.Cs20 = 3.24e-10;     % [F]  Grounded sheath S2 equivalent capacitance constant term
    
    par.Ii1  = 0.0021;       % [A]  Ion current in sheath S1
    par.Ii2  = 0.0064;       % [A]  Ion current in sheath S2
    
    par.Ie10 = 0.91;         % [A]  Electron current constant term in sheath S1
    par.Ie20 = 2.74;         % [A]  Electron current constant term in sheath S2
    
    par.Rpl  = 15.51;        % [Ohm] Plasma bulk equivalent resistance
    par.Lpl  = 2.59e-7;      % [H]   Plasma bulk equivalent inductance
    
    % ---- Geometry / circuit / operating conditions ----
    par.AE_cm2 = 100;        % [cm^2] Powered electrode area
    par.AG_cm2 = 300;        % [cm^2] Ground electrode area
    par.AE     = par.AE_cm2 * 1e-4;   % [m^2]
    par.AG     = par.AG_cm2 * 1e-4;   % [m^2]
    
    par.lbulk_cm = 5.7;      % [cm]   Plasma bulk length
    par.lbulk    = par.lbulk_cm * 1e-2;  % [m]
    
    par.Rrf   = 50;          % [Ohm] Equivalent source internal resistance
    par.Rm    = 0.5;         % [Ohm] Equivalent series resistance
    par.Rstray= 0.5;         % [Ohm] Stray parallel branch resistance
    par.Cstray= 200e-12;     % [F]   Stray parallel branch capacitance (200 pF)
    
    par.Lm2   = 1500e-9;     % [H]   Series matching inductance (1500 nH)
    % Cm1   = 1550e-12;    % [F]   Shunt matching capacitance (1550 pF)
    % Cm2   = 175e-12;     % [F]   Series matching capacitance (175 pF)
    
    par.p_mTorr = 5;         % [mTorr] Gas pressure ~= 0.667Pa
    par.p       = 0.667;
    par.Tg      = 300;       % [K]     Gas temperature
    par.Te      = 4.75;      % (as given) Electron temperature
    par.ne      = 1.25e15;   % [m^-3]  Electron density (4.5 W)
    
    par.Vrf   = 1;         % [V]   RF source voltage amplitude (25 W)
    % par.frf   = 13.56e6;     % [Hz]  RF source frequency
    % par.wrf   = 2*pi*par.frf;% [rad/s]


    %% ---------------- unpack states ----------------
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
    
    vsrc = par.Vrf*cos(wrf.*t);

    dQ1 = (-1./Cm1 .* Q1 + 1./Cm1 .* Q2 - vsrc) ./ par.Rrf;
    dQ2 = I2;
    dQ3 = I3;
    dI2 = (-(par.Rm+par.Rstray) .* dQ2 + par.Rstray .* dQ3 + 1./Cm1 .* Q1 ...
          -(1./par.Cstray + 1./Cm1 + 1./Cm2) .* Q2 + 1./par.Cstray .* Q3) ./ par.Lm2;
    dI3 = (par.Rstray .* dQ2 - (par.Rstray+par.Rpl) .* dQ3 + 1./par.Cstray .* Q2 ...
          -1./par.Cstray .* Q3 - Qs1.^2./(2*par.Cs10^2*par.Te) + Qs2.^2./(2*par.Cs20^2*par.Te)) ./ par.Lpl;
    dQs1 = dQ3  - par.Ii1 + par.Ie10 .* exp(-(Qs1./(sqrt(2)*par.Cs10*par.Te)).^2);
    dQs2 = -dQ3 - par.Ii2 + par.Ie20 .* exp(-(Qs2./(sqrt(2)*par.Cs20*par.Te)).^2);
    
    % ---------------- pack ----------------
    dx = [dQ1; dQ2; dQ3; dI2; dI3; dQs1; dQs2];

end
