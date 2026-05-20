function X_MF = solve_HOrder(X,f,rho,k,Np)
d = length(rho);
num = length(f);
h1 = -k(1):1/Np(1):k(1);
H = h1;
for ii = 2:d
    hii = -k(ii):1/Np(ii):k(ii);
    H = [kron(H,ones(1,length(hii)));kron(ones(1,length(H(1,:))),hii)];
end
f_all = rho'*H;

X_MF = [];
for ii = 1:num
    f_temp = abs(f_all - f(ii));
    [min_val, idx] = min(f_temp); 
    if f(ii)~=0 && min_val/f(ii)>1e-1
    else
       X_MF = [X_MF;H(:,idx)',X(ii,:)];
    end
end
end
