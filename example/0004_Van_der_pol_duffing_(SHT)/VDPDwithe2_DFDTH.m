function DFDTH = VDPDwithe2_DFDTH(x,dx,p,theta)
e1      = 2;
e2      = 1;
%----------------------------------------------------------
theta1 = theta(1,:);
theta2 = theta(2,:);
%----------------------------------------------------------
DFDTH = zeros(1,2,numel(theta1));
DFDTH(1,1,:)    = e1*cos(theta1);
DFDTH(1,2,:)    = e2*cos(theta2);
end