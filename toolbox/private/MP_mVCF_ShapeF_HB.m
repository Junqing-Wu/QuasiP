function [xi1,xi2,Gam0,iGam0,wGam0,U,pS,pwS] = MP_mVCF_ShapeF_HB(k,p,S,wS)
%MP_MVCF_SHAPEF_HB  Construct the periodic harmonic-balance basis for m-VCF.
%
% The basis contains one constant function and cosine/sine pairs. Its period is
% 2*P*pi, where P is the period-multiplication number p. The matrices xi1 and
% xi2 represent first- and second-derivative operators in coefficient space;
% Gam0 and iGam0 map between U-domain coefficients and S-domain samples.
% wGam0 evaluates the same basis on the dense S_wave grid.
%
% On input:
%   k  : Maximum harmonic order or explicit harmonic vector.
%   p  : Period-multiplication number P_i.
%   S  : Coarse sample count S_i on [0,2*pi).
%   wS : Dense sample count S_wave(i).
%
% On output:
%   xi1  : First-derivative matrix in coefficient space.
%   xi2  : Second-derivative matrix xi1*xi1.
%   Gam0 : Coefficient-to-sample matrix on the coarse grid.
%   iGam0: Sample-to-coefficient matrix on the coarse grid.
%   wGam0: Coefficient-to-sample matrix on the dense grid.
%   U    : Number of harmonic-balance basis functions.
%   pS   : Period-weighted coarse sample count p*S.
%   pwS  : Period-weighted dense sample count p*wS.
%
% The sampled interval is tau in [0,2*p*pi). For a scalar harmonic order k,
% the harmonics are 1/p,2/p,...,k; for an explicit vector, k is used directly.
if length(k) == 1
    h = 1/p:1/p:k;
else
    assert(p == 1, ...
      '%s: p ~= 1, but k is a vector in HB', mfilename);
    h = k;
end
%-------------------------------------------------------------------------
L   = length(h);
U   = 2*L + 1;
pS  = p*S;
pwS = p*wS;
%-------------------------------------------------------------------------
% Basis differentiation in coefficient space.
v1 = zeros(1,U-1);
v1(2:2:end) = h;
v2 = -v1;
xi1 = diag(v1,1) + diag(v2,-1);
xi2 = xi1*xi1;
%-------------------------------------------------------------------------
% Coarse-grid basis and inverse transform.
dtaui = 2*pi/S;
taui  = (0 : dtaui : 2*p*pi-dtaui)';
Gam0 = zeros(pS,U);
Gam0(:,1) = ones(pS,1);
Gam0(:,2:2:end) = cos(taui*h);
Gam0(:,3:2:end) = sin(taui*h);
iGam0 = Gam0';
iGam0(2:end,:) = 2*iGam0(2:end,:);
iGam0 = iGam0/pS;
%-------------------------------------------------------------------------
% Dense-grid basis for post-processing and nonlinear-force evaluation.
wdtaui = 2*pi/wS;
wtaui  = (0 : wdtaui : 2*p*pi-wdtaui)';
wGam0 = zeros(pwS,U);
wGam0(:,1) = ones(pwS,1);
wGam0(:,2:2:end) = cos(wtaui*h);
wGam0(:,3:2:end) = sin(wtaui*h);
end