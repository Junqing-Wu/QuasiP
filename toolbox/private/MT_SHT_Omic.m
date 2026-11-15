function [Omic,dOmic_dOm] = MT_SHT_Omic(PHI,Om1,Rho,d,Omt,Omnot)
%MT_SHT_OMIC  Construct the Fourier phase-shift matrix for SHT closure.
%
% The multi-frequency shooting closure requires a coefficient transformation
% for the shift theta -> theta-2*pi*Rho. With kQ=PHI.kQ and
% P=cell2mat(PHI.P), the phase angles used by the code are
%   lambda = 2*P(Omt)*pi*Rho'*kQ.
% Omic is assembled so that
%   Phi(theta-2*pi*Rho) = Phi(theta)*Omic.
% The cell array dOmic_dOm contains the derivatives needed in the Jacobian.
%
% On input:
%   PHI   : Shooting discretization data, including kQ and period factors P.
%   Om1   : Selected integration frequency Om(Omt).
%   Rho   : Ratios Om(Omnot)/Om1.
%   d     : Number of base frequencies.
%   Omt   : Index of the integration frequency.
%   Omnot : Indices of the remaining frequencies.
%
% On output:
%   Omic       : Sparse Fourier phase-shift matrix.
%   dOmic_dOm  : Derivatives of Omic with respect to all entries of Om.
kQ = PHI.kQ;
P  = cell2mat(PHI.P);
lamda = 2*P(Omt)*pi*Rho'*kQ;

v0   = [1,kron(cos(lamda),[1,1])];
v1   = kron(sin(lamda),[0,1]);
Omic = sparse(diag(v0)+diag(-v1,1)+diag(v1,-1));

dOmic_dOm = cell(1,d);
for  ii = 1:d-1
     v0 = [0,zeros(1,2*size(kQ,2))];
     v1 = kron(kQ(ii,:),[0,1]);
     D = diag(v0)+diag(-v1,1)+diag(v1,-1);
     dOmic_dOm{Omnot(ii)} = 2*P(Omt)*pi/Om1*D*Omic;
end
v1 = kron(lamda,[0,1]);
D = diag(v0)+diag(-v1,1)+diag(v1,-1);
dOmic_dOm{Omt} = -1/Om1*D*Omic;
end