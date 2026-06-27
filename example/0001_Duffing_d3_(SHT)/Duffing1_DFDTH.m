function DFDTH = Duffing1_DFDTH(x,p,theta)
%----------------------------------------------------------
theta1 = theta(1,:);
theta2 = theta(2,:);
theta3 = theta(3,:);
%----------------------------------------------------------
% zeta = p(1,:);
% gama = p(2,:);
f1   = p(3,:);
f2   = p(4,:);
%----------------------------------------------------------
DFDTH = zeros(1,3,numel(theta1));
DFDTH(1,1,:) = f1.*cos(theta1);
DFDTH(1,2,:) = f2.*cos(theta2);
DFDTH(1,3,:) = 5*cos(theta3);
end