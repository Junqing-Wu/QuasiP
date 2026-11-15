function PHI = MT_SHT_COM_PHI_Supplement(data)
% Supplementary SHT relation:
% The supplement adds the phase-shift and sensitivity data needed to evaluate
% R_M after the last trajectory segment.
%MT_SHT_COM_PHI_SUPPLEMENT  Complete PHI data for MT-SHT residual evaluation.
%
% Add model- and method-dependent quantities that cannot be constructed by
% the generic MT-SHT mesh builder alone, including trajectory sensitivities,
% solver callbacks, phase-direction data, and period-shift information.
%
% On input:
%   data : QuasiP data containing the model callbacks and initial PHI fields.
%
% On output:
%   PHI  : Updated PHI structure ready for MT_SHT_QuasiP.
PHI      = data.PHI;
PHI.N    = data.ODE*PHI.n;
PHI.nL   = PHI.N*PHI.L;
%-------------------------------------------------------------
local_HT = cell(1,PHI.M);
local_MF = cell(1,PHI.M);
for jj = 1:PHI.M
    local_HT{jj} = (jj-1)*PHI.pS2d+1:jj*PHI.pS2d;
    local_MF{jj} = (jj-1)*PHI.U+1:jj*PHI.U;
end
PHI.local_HT = local_HT;
PHI.local_MF = local_MF;
%-------------------------------------------------------------
Local_HT = cell(1,PHI.M);
Local_MF = cell(1,PHI.M);
for jj = 1:PHI.M
    Local_HT{jj} = [];
    Local_MF{jj} = [];
    for nn = 1:PHI.N
        Local_HT{jj} = [Local_HT{jj},(nn-1)*PHI.L_HT+local_HT{jj}];
        Local_MF{jj} = [Local_MF{jj},(nn-1)*PHI.L_MF+local_MF{jj}];
    end
end
PHI.Local_HT = Local_HT;
PHI.Local_MF = Local_MF;
if data.ode_opts.d > 1
    PHI.EiG      = sparse( kron(eye(PHI.N),PHI.iGam) );
    PHI.EG       = sparse( kron(eye(PHI.N),PHI.Gam ) );
end
if data.ODE == 1 && ~PHI.MATLABsolver
    PHI.ODEsolver = @ode45;
    PHI.MATLABsolver = true;
    PHI.IsFixedStep  = false;
    disp('First-order systems can only use MATLABsolver! >>>ODE45 is used by default<<<');
end
end