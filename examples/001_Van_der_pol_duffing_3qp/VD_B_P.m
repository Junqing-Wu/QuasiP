function B_p = VD_B_P(p)

mu  = p(1);
Om0 = 2;
alp = 0.5;
eps1 = 2;
eps2 = 1;
%%
B_p = zeros(2,2,numel(p));

B_p(:,:,1) = [-1  0;
               0  0]; 
end