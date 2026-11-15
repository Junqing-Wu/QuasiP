function [x_mf,A] = HT2MF(x_ht,pwS,pwS_prod)
%HT2MF  Convert hyper-time samples to multi-frequency coefficients.
%
% Reshape x_ht on the fine grid pwS, apply the multidimensional FFT, normalize
% by pwS_prod, and sort coefficient energies. A is the retained coefficient
% index list used when constructing an initial Xk_MF representation.
%
% On input:
%   x_ht     : Hyper-time sample array.
%   pwS      : Number of fine samples in each hyper-time direction.
%   pwS_prod : Product of the entries of pwS.
%
% On output:
%   x_mf     : Flattened multi-frequency coefficient array.
%   A        : Indices sorted by energy and filtered by the relative threshold.

x_ht = reshape(x_ht,pwS);
X_mf = fftn(x_ht)/pwS_prod;
x_mf = reshape(X_mf,1,pwS_prod);
[energy, A] = sort(abs(x_mf).^2, 'descend'); 
energy = energy/energy(1);
idx1 = energy >= 1e-14;
A = A(idx1);
end