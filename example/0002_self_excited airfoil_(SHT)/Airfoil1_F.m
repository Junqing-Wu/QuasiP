function F = Airfoil1_F(x,p,theta)
kbeta3 = 10;
q3 = x(3,:);
S = length(q3);
%----------------------------------------------------------
F = zeros(3,S);
F(3,:) = - kbeta3*q3.^3;
end
