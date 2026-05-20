function [xi1,xi2,Gam0,iGam0,wGam0,U,pS,pwS] = MP_mVCF_ShapeF_HB(k,p,S,wS)
if length(k) == 1
    h = 1/p:1/p:k;
else
    assert(p == 1, ...
    '%s: p ~= 1, but k is a vector in HB');
    h = k;
end
%----------------------------------------------------------
L   = length(h);
U   = 2*L + 1;
pS  = p*S;
pwS = p*wS;
%----------------------------------------------------------
v1 = zeros(1,U-1);          
v1(2:2:end) = h;            
v2 = -v1; 
xi1 = diag(v1,1) + diag(v2,-1);
xi2 = xi1*xi1;
%----------------------------------------------------------
dtaui = 2*pi/S;
Taui  = 2*p*pi;
taui  = (0 : dtaui : Taui-dtaui)';
Gam0 = zeros(pS,U);
Gam0(:,1) = ones(pS,1);   
Gam0(:,2:2:end) = cos(taui*h);            
Gam0(:,3:2:end) = sin(taui*h);
iGam0 = Gam0';
iGam0(2:end,:) = 2*iGam0(2:end,:);             
iGam0 = iGam0/(pS);                           
%----------------------------------------------------------
wdtaui = 2*pi/wS;
Taui   = 2*p*pi;
wtaui  = (0 : wdtaui : Taui-wdtaui)';
wGam0 = zeros(pwS,U);
wGam0(:,1) = ones(pwS,1);   
wGam0(:,2:2:end) = cos(wtaui*h);            
wGam0(:,3:2:end) = sin(wtaui*h);
end

% y = cos(0.5*taui) + sin(3*taui) + cos(taui);
% dy = -0.5*sin(0.5*taui) + 3*cos(3*taui) - sin(taui);
% ddy = -0.5*0.5*cos(0.5*taui) - 9*sin(3*taui) - cos(taui);
% 
% Y = iGam0*y;
% yy = Gam0*Y;
% wyy = wGam0*Y;
% dyy = Gam0*xi1*Y;
% wdyy = wGam0*xi1*Y;
% ddyy = Gam0*xi2*Y;
% wddyy = wGam0*xi2*Y;
% 
% figure(1)
% plot(taui, y, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, yy, 'r--', 'LineWidth', 1.5);
% plot(wtaui, wyy, 'b--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','wt');

% figure(2)
% plot(taui, dy, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, dyy, 'r--', 'LineWidth', 1.5);
% plot(wtaui, wdyy, 'b--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','wt');

% figure(3)
% plot(taui, ddy, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, ddyy, 'r--', 'LineWidth', 1.5);
% plot(wtaui, wddyy, 'b--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','wt');
