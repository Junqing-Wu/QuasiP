function [data, y, J] = MT_SHT_QuasiP(prob, data, u)
%MT_SHT_QUASIP  Integrate the trajectory family and assemble FSE-SHT residuals.
%
% This is the main multi-trajectory shooting residual. It reshapes the COCO
% unknown vector into initial trajectory data, integrates every segment, stores
% terminal states and sensitivities in PSI, and then calls the hyper-time or
% multi-frequency boundary routine selected by PHI.type.
%
% On input:
%   data : QuasiP model callbacks, PHI shooting settings, and post-processing options.
%   u    : COCO unknown vector. Its first PHI.nL entries contain the shooting
%          coefficients, followed by d base frequencies and model parameters.
%
% On output:
%   data : Updated data; amplitude and stability trajectory data may be stored.
%   y    : FSE-SHT boundary residual returned to COCO.
%   J    : Jacobian of y with respect to the active continuation variables.
%
% The trajectory equation is integrated along phi1=Om(Omt)*t. The remaining
% frequencies enter through Rho=Om(Omnot)/Om1. The selected PHI.P period factors
% determine the phase shift used when closing the torus.
PHI    = data.PHI;                  % PHI basis, transform, and differentiation data
d      = data.ode_opts.d;           % Number d of response base frequencies
MATLABsolver = PHI.MATLABsolver;    % 1：y' = 1/Om1*g(y,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)  
                                    % 0: M*a + D*u + K*q = f(q,Om1*u,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
                                    %    M = Om1^2*A; D = Om1*B 
%----------------------------------------------------------
N      = PHI.N;                     % First-order state dimension N
NL     = PHI.nL;                    % Total coefficient count nL=N*L
M_SHT  = PHI.M;                 % Number M of shooting segments
a      = data.pdim;                 % Number a of model parameters p
Omt    = PHI.Omt;                   % Integration direction Omt; phi1=Om(Omt)*t and period 2*P(Omt)*pi
%----------------------------------------------------------
switch PHI.type
    case 'HT'  % SHT_HT hyper-time shooting 
               % d=1 is evaluated in HT-SHT
        ys  = reshape(u(1:NL),[],N);
    case 'MF'  % SHT_MF multi-frequency shooting
        Ys  = reshape(u(1:NL),[],N*M_SHT);
        ys  = reshape(PHI.Gam*Ys,[],N);
        Ys  = reshape(u(1:NL),[],N);        
end
Om  = u(1+NL:d+NL);
p   = u(d+NL+1:end);
%----------------------------------------------------------
phis  = PHI.Varphi;        % Initial phases of shooting segments
phi1s = PHI.Varphi1_sta;   % Segment start values in phi1
phi1e = PHI.Varphi1_end;   % Segment end values in phi1
dphi1  = PHI.dphi1;    % Fixed integration step in phi1 
%==========================================================
%>>>>>>>>>>>>>>>>>>>>>>>>> TI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%==========================================================
if PHI.isvar
    PSI = zeros(PHI.L_HT,N+N^2+N*d+N*a);
else
    PSI = zeros(PHI.L_HT,N);
end
if data.amp_opts.isamp
    if MATLABsolver
        phi_y_dy = cell(PHI.L_HT,1);
    else
        phi_y_dy = zeros(PHI.num_phi1_taj,1+N,PHI.L_HT);
    end
end
% The solver needs model and discretization data, not the preceding
% residual's trajectory history. Keep those arrays out of worker transfers.
solverData = data;
history = intersect(fieldnames(solverData),{'phi_y_dy','ye_ys'});
if ~isempty(history)
    solverData = rmfield(solverData,history);
end
if data.par_opts.parallel
    if MATLABsolver
        % y' = 1/Om1*g(y,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
        parfor ii = 1 : PHI.L_HT
            ys_ii     = ys(ii,:);
            phis_ii   = phis(ii,:);
            phi1s_ii  = phis_ii(Omt);
            idx_match = (phi1s == phi1s_ii);
            phi1e_ii  = phi1e(idx_match);
            [varphi1,Psi,dy]   = MT_SHT_MATLABsolver1(solverData,ys_ii,phis_ii,phi1s_ii,phi1e_ii,Om,p);
            PSI(ii,:) = Psi(end,:);
            if data.amp_opts.isamp
                phi_y_dy{ii,1} = [varphi1,Psi(:,1:N),dy];
            end
        end
    else
        % M*a + D*u + K*q = f(q,Om1*u,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
        parfor ii = 1 : PHI.L_HT
            ys_ii     = ys(ii,:);
            phis_ii   = phis(ii,:);
            phi1s_ii  = phis_ii(Omt);
            idx_match = (phi1s == phi1s_ii);
            phi1e_ii  = phi1e(idx_match);
            [varphi1,Psi]   = PHI.ODEsolver(solverData,ys_ii(:),phis_ii(:),phi1s_ii,phi1e_ii,dphi1,Om,p);
            PSI(ii,:) = Psi(end,:);
            if data.amp_opts.isamp
                phi_y_dy(:,:,ii) = [varphi1,Psi(:,1:N)];
            end
        end
    end
else
    if MATLABsolver
        % y' = 1/Om1*g(y,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
        for ii = 1 : PHI.L_HT
            ys_ii     = ys(ii,:);
            phis_ii   = phis(ii,:);
            phi1s_ii  = phis_ii(Omt);
            idx_match = (phi1s == phi1s_ii);
            phi1e_ii  = phi1e(idx_match);
            [varphi1,Psi,dy]   = MT_SHT_MATLABsolver1(solverData,ys_ii,phis_ii,phi1s_ii,phi1e_ii,Om,p);
            PSI(ii,:) = Psi(end,:);
            if data.amp_opts.isamp
                phi_y_dy{ii,1} = [varphi1,Psi(:,1:N),dy];
            end
        end
    else
        % M*a + D*u + K*q = f(q,Om1*u,p,phi1,phi2+rho2*phi1,...,phie+rhoe*phi1)
        for ii = 1 : PHI.L_HT
            ys_ii     = ys(ii,:);
            phis_ii   = phis(ii,:);
            phi1s_ii  = phis_ii(Omt);
            idx_match = (phi1s == phi1s_ii);
            phi1e_ii  = phi1e(idx_match);
            [varphi1,Psi]   = PHI.ODEsolver(solverData,ys_ii(:),phis_ii(:),phi1s_ii,phi1e_ii,dphi1,Om,p);
            PSI(ii,:) = Psi(end,:);
            if data.amp_opts.isamp
                phi_y_dy(:,:,ii) = [varphi1,Psi(:,1:N)];
            end
        end
    end
end
%==========================================================
%>>>>>>>>>>>>>>>>>>>>>>>>> BC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%==========================================================
switch PHI.type
    case 'HT'
        [y,J] = MT_SHT_HT_QuasiP(data,PSI,ys,Om);
    case 'MF'
        [y,J] = MT_SHT_MF_QuasiP(data,PSI,Ys,Om);
end

if data.amp_opts.isamp && max(abs(y))<1e-5 
        data.phi_y_dy  = phi_y_dy;
end
if data.sta_opts.issta && PHI.isvar && max(abs(y))<1e-5
    data.ye_ys = PSI(:,N+1:N+N^2);
end
end
