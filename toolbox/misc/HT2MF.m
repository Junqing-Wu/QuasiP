function [x_mf,A] = HT2MF(x_ht,pwS,pwS_prod)
x_ht = reshape(x_ht,pwS);
X_mf = fftn(x_ht)/pwS_prod;
x_mf = reshape(X_mf,1,pwS_prod);
[energy, A] = sort(abs(x_mf).^2, 'descend'); 
energy = energy/energy(1);
idx1 = energy >= 1e-14;
A = A(idx1);
end