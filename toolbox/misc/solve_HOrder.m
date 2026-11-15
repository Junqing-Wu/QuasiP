function Xk_MF = solve_HOrder(Xc,f,rho,H,P)
%SOLVE_HORDER  Construct an Xk_MF initial guess from spectral peaks.
%
% Build candidate harmonic combinations k_tilde from the base-frequency ratios
% rho and the requested harmonic orders H/P. Match rho''*k_tilde to the
% detected frequencies f and retain the corresponding coefficient rows Xc.
%
% On input:
%   Xc   : Spectral coefficient rows associated with f.
%   f    : Detected physical frequencies.
%   rho  : Base-frequency ratios.
%   H    : Maximum harmonic order in each frequency direction.
%   P    : Period-multiplication numbers used in the harmonic grid.
%
% On output:
%   Xk_MF.k_wave : Selected harmonic combinations k_tilde.
%   Xk_MF.Xc     : Fourier coefficient rows retained for those combinations.

d = length(rho);
num = length(f);
h1 = -H(1):1/P(1):H(1);
kp = h1;
for ii = 2:d
    hii = -H(ii):1/P(ii):H(ii);
    kp = [kron(kp,ones(1,length(hii)));kron(ones(1,length(kp(1,:))),hii)];
end
f_all = rho'*kp;

Xk_MF.k_wave = [];
Xk_MF.Xc = [];
for ii = 1:num
    f_temp = abs(f_all - f(ii));
    [min_val, idx] = min(f_temp); 
    if f(ii)~=0 && min_val/f(ii)>1e-1
    else
       Xk_MF.k_wave = [Xk_MF.k_wave; kp(:,idx)'];
       Xk_MF.Xc = [Xk_MF.Xc; Xc(ii,:)];
    end
end
end