function [x_mf,A] = HT2MF(x_ht,pwS,pwS_prod)
x_ht = reshape(x_ht,pwS);
X_mf = fftn(x_ht)/pwS_prod;
x_mf = reshape(X_mf,1,pwS_prod);

A = find( abs(x_mf) >= max(abs(x_mf))*10^(-8) );
end