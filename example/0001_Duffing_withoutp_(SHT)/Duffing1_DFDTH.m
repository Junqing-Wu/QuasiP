function DFDTH = Duffing1_DFDTH(x,p,theta)
%----------------------------------------------------------
theta1 = theta(1,:);
theta2 = theta(2,:);
%----------------------------------------------------------
% zeta = 0.1;
% gama = 0.2;
f1   = 5;
f2   = 5;
%----------------------------------------------------------
DFDTH = zeros(1,2,numel(theta1));
DFDTH(1,1,:) = f1.*cos(theta1);
DFDTH(1,2,:) = f2.*cos(theta2);
end