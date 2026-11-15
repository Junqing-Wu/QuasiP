function C = MT_SHT_INITIAL_C_byPHI(data,Xk_MF,Om)
%MT_SHT_INITIAL_C_BYPHI  Build initial shooting coefficients from Xk_MF.
%
% Convert the user multi-frequency representation to the initial trajectory
% or Fourier coefficient ordering required by the selected MT-SHT type. OM is
% used to construct the phase-aligned initial sections when required.
%
% On input:
%   data : QuasiP data, initial multi-frequency representation, and
%   Xk_MF : QuasiP data, initial multi-frequency representation, and
%   Om : QuasiP data, initial multi-frequency representation, and
%                   base frequencies.
%
% On output:
%   C : Initial MT-SHT coefficient vector.
PHI     = data.PHI;
d       = data.ode_opts.d;
%--------------------------------------------------------------------------------------------------------------------
Omt     = PHI.Omt;
Omnot   = PHI.Omnot;
%--------------------------------------------------------------------------------------------------------------------
Rho     = Om(Omnot)./Om(Omt);
Varphi  = PHI.Varphi;
if d == 1
    Theta   = Varphi;
else
    Rhophi1 = Varphi(:,Omt).*Rho';
    Theta   = Varphi;
    Theta(:,Omnot) = Theta(:,Omnot) + Rhophi1;
end
%--------------------------------------------------------------------------------------------------------------------
k_h = Xk_MF.k_wave';
Q   = Xk_MF.Xc;
%--------------------------------------------------------------------------------------------------------------------
H_IDF = exp(1i*Theta*k_h);
kOm   = 1i*diag(Om'*k_h);
%--------------------------------------------------------------------------------------------------------------------
qs    = real(H_IDF*Q);
vs    = real(H_IDF*kOm*Q);
us    = real(H_IDF*kOm*Q)/Om(Omt);
%--------------------------------------------------------------------------------------------------------------------
if d == 1
    if data.ODE == 1
        C = reshape(qs,[],1);
    else
        if PHI.MATLABsolver
            C = [reshape(qs,[],1);reshape(vs,[],1)];
        else
            C = [reshape(qs,[],1);reshape(us,[],1)];
        end
    end
else
    pS2d = PHI.pS2d;
    qs   = reshape(qs,pS2d,[]);
    vs   = reshape(vs,pS2d,[]);
    us   = reshape(us,pS2d,[]);
    Qs   = PHI.iGam*qs;
    Vs   = PHI.iGam*vs;
    Us   = PHI.iGam*us;
    switch PHI.type
        case 'HT'
            if data.ODE == 1
                C = reshape(qs,[],1);
            else
                if PHI.MATLABsolver
                    C = [reshape(qs,[],1);reshape(vs,[],1)];
                else
                    C = [reshape(qs,[],1);reshape(us,[],1)];
                end
            end
        case 'MF'
            if data.ODE == 1
                C = reshape(Qs,[],1);
            else
                if PHI.MATLABsolver
                    C = [reshape(Qs,[],1);reshape(Vs,[],1)];
                else
                    C = [reshape(Qs,[],1);reshape(Us,[],1)];
                end
            end
    end
end
end