function PHI = MT_SHT_COM_PHI(PHI_INITIAL)
% SHT discretization relation:
%   Y_s -> y_s = [I_N kron Gam]*Y_s,
%   y_e -> Y_e = [I_N kron iGam]*y_e.
% PHI stores the segment maps required by the MF and HT boundary residuals.
%MT_SHT_COM_PHI  Build the multi-trajectory shooting discretization.
%
% Construct the trajectory mesh, Fourier index combinations, sample/coefficient
% transforms, segment location maps, and period-dependent phase-shift data used
% by the MT-SHT residual. The period factors PHI.P determine the length and
% phase relation of each hyper-time direction.
%
% On input:
%   PHI_INITIAL : MT-SHT settings including d, n, method type, harmonic orders,
%                 period factors, mesh counts, and the selected shooting direction.
%
% On output:
%   PHI : Complete MT-SHT discretization structure.
PHI = PHI_INITIAL;
d   = PHI.d;
%---------------------------------------------------------------------
defaults.ODEsolver      = @ode45;
defaults.Omt            = 1;
defaults.M              = 1;
defaults.L1             = 2^8;
defaults.IsFixedStep    = false;
defaults.isvar          = true;
defaults.ODEsolver_opt  = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-8, 'NormControl', 'on');
defaults.ODE2slover_opt = struct('ItMX',10,'RelTol', 1.0e-6,'alpha',0,'rhoinf',0.9);
def_fields = fieldnames(defaults);
for f = def_fields'
    fname = f{1};
    if ~isfield(PHI, fname)
        PHI.(fname) = defaults.(fname);
    end
end
%---------------------------------------------------------------------
ODEsolver = PHI.ODEsolver;
matlab_solver_list = {@ode45, @ode23, @ode113, @ode15s, @ode23s, @ode23t, @ode23tb};
PHI.MATLABsolver = any(cellfun(@(solver) isequal(ODEsolver, solver), matlab_solver_list));
if PHI.MATLABsolver
    PHI.IsFixedStep = false;
else
    PHI.IsFixedStep = true;
end
%---------------------------------------------------------------------
PHI.Omnot      = [PHI.Omt+1:d,1:PHI.Omt-1];
knot       = PHI.k(PHI.Omnot);
Pnot       = PHI.P(PHI.Omnot);
Snot       = PHI.S(PHI.Omnot);
S_wavenot  = PHI.S_wave(PHI.Omnot);
if d == 1
    PHI.kQ       = [];    
    PHI.U        = 1;
    PHI.pS2d     = 1;
    PHI.Varphi2d = [];
    PHI.type     = 'HT';   % d=1 uses HT-SHT by default
else
    %---------------------------------------------------------------------
    PHI.kQ                     = Com_kQ(knot,Pnot,d);
    PHI.U                      = 2*length(PHI.kQ(1,:))+1;
    [PHI.Varphi2d,PHI.pSnot]   = Com_phi(Snot,Pnot);
    PHI.pS2d                   = prod(cell2mat(PHI.pSnot));

    [PHI.Gam,PHI.iGam]         = Com_Gam(PHI.kQ,PHI.U,PHI.Varphi2d,cell2mat(PHI.pSnot));
    PHI.iGTiG                  = transpose(PHI.iGam)*PHI.iGam;

    [Varphi_wave2d,PHI.pwSnot] = Com_phi(S_wavenot,Pnot);
    [PHI.wGam,PHI.iwGam]       = Com_Gam(PHI.kQ,PHI.U,Varphi_wave2d,cell2mat(PHI.pwSnot));
    PHI.wGamiGam               = PHI.wGam*PHI.iGam;
    %---------------------------------------------------------------------
    v1  = zeros(1,PHI.U-1);
    v1(2:2:end) = 1;
    v2     = -v1;
    Xi1    = diag(v1,1)+diag(v2,-1);
    PHI.XI1_i    = cell(1,d);
    PHI.iGTXI1iG = cell(1,d);
    PHI.OkQ      = [zeros(d-1,1),kron(PHI.kQ,[1,1])];
    for ii = 1:d-1
        PHI.XI1_i{PHI.Omnot(ii)}    = sparse( Xi1.*PHI.OkQ(ii,:) );
        PHI.iGTXI1iG{PHI.Omnot(ii)} = transpose(PHI.iGam)*PHI.XI1_i{PHI.Omnot(ii)}*PHI.iGam;
    end
end
%---------------------------------------------------------------------
if PHI.IsFixedStep
    if PHI.S_wave{PHI.Omt} >= PHI.L1
        n_node    = PHI.P{PHI.Omt}*PHI.L1 + 1;
        mm        = PHI.S_wave{PHI.Omt}/PHI.L1;
        n_interp  = PHI.P{PHI.Omt}*PHI.S_wave{PHI.Omt} + 1;
        [A,B,C,D] = cubic_spline(PHI.P{PHI.Omt}*PHI.L1, PHI.P{PHI.Omt});
        SPLine    = zeros(n_interp, n_node);
        SPLine(1:mm:end, :) = eye(n_node);
        PHI.dphi1           = 2*pi/PHI.L1;
        for seg = 1:PHI.P{PHI.Omt}*PHI.L1
            for ii = 1:mm-1
                deltexi = ii/mm * PHI.dphi1;
                idx = (seg-1)*mm + ii + 1;
                SPLine(idx, :) = A(seg, :) + B(seg, :)*deltexi + C(seg, :)*deltexi^2 + D(seg, :)*deltexi^3;
            end
        end
        PHI.SPLine = sparse(SPLine);
    else
        mm = PHI.L1/PHI.S_wave{PHI.Omt};
        SPLine = eye(PHI.P{PHI.Omt}*PHI.L1+1);
        PHI.SPLine = sparse(SPLine(1:mm:end,:));
    end
end
PHI.Varphi1      = 0:2*pi/PHI.L1:2*PHI.P{PHI.Omt}*pi;
PHI.wVarphi1      = 0:2*pi/PHI.S_wave{PHI.Omt}:2*PHI.P{PHI.Omt}*pi;
%---------------------------------------------------------------------
Dphi1             = 2*PHI.P{PHI.Omt}*pi/PHI.M;
Varphi1_se        = 0:Dphi1:2*PHI.P{PHI.Omt}*pi;
PHI.Varphi1_sta   = Varphi1_se(1:end-1);
PHI.Varphi1_end   = Varphi1_se(2:end);
%---------------------------------------------------------------------
Varphi              = zeros(PHI.M*PHI.pS2d,d);
Varphi(:,PHI.Omt)   = kron(transpose(PHI.Varphi1_sta),ones(PHI.pS2d,1));
Varphi(:,PHI.Omnot) = kron(ones(PHI.M,1),PHI.Varphi2d);
PHI.Varphi          = Varphi;             % Initial trajectory phases
PHI.L_HT            = PHI.M*PHI.pS2d; % Unknown count for one state in HT-SHT
PHI.L_MF            = PHI.M*PHI.U;    % Unknown count for one state in MF-SHT
switch PHI.type
    case 'HT'
        PHI.L = PHI.L_HT;
    case 'MF'
        PHI.L = PHI.L_MF;
end
PHI.dphi1  = 2*pi/PHI.L1;
PHI.num_phi1_taj = PHI.P{PHI.Omt}*PHI.L1/PHI.M + 1;
%---------------------------------------------------------------------
OmOd     = [PHI.Omt,PHI.Omnot];
[~, pos] = sort(OmOd);
P_omod   = cell(1,d);
wS_omod  = cell(1,d);
pwS_omod = cell(1,d);
for ii = 1:d
    P_omod{ii}   = PHI.P{OmOd(ii)};
    wS_omod{ii}  = PHI.S_wave{OmOd(ii)};
    pwS_omod{ii} = PHI.P{OmOd(ii)}*PHI.S_wave{OmOd(ii)};
end
PHI.OmOd     = OmOd;
PHI.pos      = pos;
PHI.P_omod   = P_omod;
PHI.wS_omod  = wS_omod;
PHI.pwS_omod = pwS_omod;
%---------------------------------------------------------------------
wS  = cell2mat(PHI.wS_omod) ;
pwS = cell2mat(PHI.pwS_omod);
p   = cell2mat(PHI.P_omod)  ;
wdtau1   = 2*pi/wS(1);
wTau1    = 2*p(1)*pi;
wtau_1im1 = (0:wdtau1:wTau1-wdtau1)';
pwS_1im1 = pwS(1);
for ii = 2:d
    wdtaui = 2*pi/wS(ii);
    wTaui  = 2*p(ii)*pi;
    wtau_i = (0:wdtaui:wTaui-wdtaui)';  
    wtau_1im1 = [ kron( ones(pwS(ii),1),wtau_1im1 ) , kron( wtau_i,ones(pwS_1im1,1) ) ];
    pwS_1im1  = pwS_1im1*pwS(ii);  
end
wtau     = wtau_1im1';
wtau     = wtau(pos,:);
PHI.wtau = wtau;
end
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
function kQ = Com_kQ(knot,Pnot,d)
h = cell(1,d-1);
N = cell(1,d-1);
for ii = 1:d-1
    k = knot{ii};
    p = Pnot{ii};
    if length(k) == 1
        h{ii} = 1/p:1/p:k;
    else
        assert(p == 1, ...
            '%s: p ~= 1, but k is a vector in HB');
        h{ii} = k;
    end
    N{ii} = length(h{ii});
end
kQ = [0,h{1}];
for ii = 2:d-1
    kQ = kQ(:,2:end);
    NQ = size(kQ,2);
    OO = zeros(ii-1,1);
    kQ  = [ [kron(ones(1,N{ii}+1),[OO,kQ]);kron([0,h{ii}],ones(1,NQ+1))],...
        [kron(ones(1,N{ii})  ,kQ)     ;kron(-h{ii}   ,ones(1,NQ))] ];
end
kQ = kQ(:,2:end);
end
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
function [Varphi,pS] = Com_phi(Snot,Pnot)
pS  = cell(1,numel(Snot));
S   = cell2mat (Snot);
P   = cell2mat (Pnot);
pS{1} = P(1)*S(1);

Varphi1IM1 = (0:2*pi/S(1):2*P(1)*pi-2*pi/S(1))';
pS_1im1 = pS{1};
for ii = 2:numel(Snot)
    pS{ii} = P(ii)*S(ii);
    Varphi_i = (0:2*pi/S(ii):2*P(ii)*pi-2*pi/S(ii))';
    Varphi1IM1 = [ kron( ones(pS{ii},1),Varphi1IM1 ) , kron( Varphi_i,ones(pS_1im1,1) ) ];
    pS_1im1 = pS_1im1*pS{ii};
end
Varphi = Varphi1IM1;
end
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
function [Gam,iGam] = Com_Gam(kQ,UQ,Varphi,pS_2d)
Gam      = zeros( prod(pS_2d), UQ );
Gam(:,1) = ones( prod(pS_2d),1 );
Gam(:,2:2:end) = cos(Varphi*kQ);
Gam(:,3:2:end) = sin(Varphi*kQ);
iGam = Gam';
iGam(2:end,:) = 2*iGam(2:end,:);
iGam = iGam/( prod(pS_2d) );
end