function B_p = VD_B_P(p)

mu  = p(1);
Om0 = p(2);
alp = p(3);
eps1 = p(4);
eps2 = p(5);
%%
B_p = zeros(2,2,numel(p));

B_p(:,:,1) = [-1  0;
               0  0]; 
end