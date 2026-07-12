function DFDTH = VDPDwithe3_DFDTH(x,dx,p,theta)
e1      = 2;
e2      = 1;
e3      = 0.5;
%----------------------------------------------------------
theta1 = theta(1,:);
theta2 = theta(2,:);
theta3 = theta(3,:);
%----------------------------------------------------------
DFDTH = zeros(1,3,numel(theta1));
DFDTH(1,1,:)    = e1*cos(theta1);
DFDTH(1,2,:)    = e2*cos(theta2);
DFDTH(1,3,:)    = e3*cos(theta3);
end