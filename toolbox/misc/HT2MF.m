function [x_mf,A] = HT2MF(x_ht,pwS,pwS_prod)
x_ht = reshape(x_ht,pwS);
X_mf = fftn(x_ht)/pwS_prod;
x_mf = reshape(X_mf,1,pwS_prod);
[a_xmf, A] = sort(abs(x_mf), 'descend'); 
sum_a = cumsum(a_xmf);
sum_a = sum_a/sum_a(end);
N_9992  = find(sum_a < 0.9992, 1, 'last');
N_99992 = find(sum_a < 0.99992, 1, 'last');

if N_9992/N_99992 < 0.25
    A = A(1:N_9992);
else
    A = A(1:N_99992);
end

end