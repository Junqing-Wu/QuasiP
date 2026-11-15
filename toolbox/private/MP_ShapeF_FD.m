function [xi1,xi2,wGam0,wGam1,U,pS,pwS] = MP_ShapeF_FD(w,p,S,wS)
%PURPOSE  Construct periodic finite-difference operators and interpolation.
%
% Build periodic first- and second-derivative matrices from the requested
% finite-difference stencil and construct interpolation operators on the dense
% S_wave grid.
%
% On input:
%   w  : Finite-difference stencil positions.
%   p  : Period-multiplication number P_i.
%   S  : Coarse sample count S_i.
%   wS : Dense sample count S_wave(i).
%
% On output:
%   xi1 : Periodic first-derivative matrix.
%   xi2 : Periodic second-derivative matrix.
%   wGam0: Dense-grid value interpolation matrix.
%   wGam1: Dense-grid first-derivative interpolation matrix.
%   U   : Number of finite-difference coefficients.
%   pS  : Period-weighted coarse sample count P_i*S_i.
%   pwS : Period-weighted dense sample count P_i*S_wave(i).
tol = 1e-6;
if abs(mod(wS, S)) > tol
    warning('S_wave is not an integer multiple of S!');
end
pS  = p*S;
U   = p*S;
pwS = p*wS;
%----------------------------------------------------------
dtaui = 2*pi/S;
%----------------------------------------------------------
m = length(w); 
TY1 = zeros(m,1);TY1(2) = 1;    
TY2 = zeros(m,1);TY2(3) = 2*1;
W = zeros(m,m);
for ii = 0:m-1
    W(ii+1,:) = w.^ii;
end
alp1 = W\TY1; alp2 = W\TY2;
%----------------------------------------------------------
for ii = 1:m
    if abs(alp1(ii))<10e-8 
        alp1(ii) = 0; 
    end
    if abs(alp2(ii))<10e-8
        alp2(ii) = 0; 
    end
end
%----------------------------------------------------------
[~,I_0] = find(w==0);
[~,I_p] = find(w >0);  
[~,I_n] = find(w <0);
xi1 = diag(ones(U,1)*alp1(I_0),0);
xi2 = diag(ones(U,1)*alp2(I_0),0);
for ii = 1:length(I_p)
    xi1 = xi1 + diag(ones(U-abs(w(I_p(ii))),1)*alp1(I_p(ii)),w(I_p(ii)));
    xi1 = xi1 + diag(ones(U-(U-abs(w(I_p(ii)))),1)*alp1(I_p(ii)),-U+w(I_p(ii)));
    xi2 = xi2 + diag(ones(U-abs(w(I_p(ii))),1)*alp2(I_p(ii)),w(I_p(ii)));
    xi2 = xi2 + diag(ones(U-(U-abs(w(I_p(ii)))),1)*alp2(I_p(ii)),-U+w(I_p(ii)));
end
for ii = 1:length(I_n)
    xi1 = xi1 + diag(ones(U-abs(w(I_n(ii))),1)*alp1(I_n(ii)),w(I_n(ii)));
    xi1 = xi1 + diag(ones(U-(U-abs(w(I_n(ii)))),1)*alp1(I_n(ii)),U+w(I_n(ii)));
    xi2 = xi2 + diag(ones(U-abs(w(I_n(ii))),1)*alp2(I_n(ii)),w(I_n(ii)));
    xi2 = xi2 + diag(ones(U-(U-abs(w(I_n(ii)))),1)*alp2(I_n(ii)),U+w(I_n(ii)));
end
xi1 = xi1/dtaui;
xi2 = xi2/dtaui^2;
%----------------------------------------------------------
[A,B,C,D] = periodic_spline(U,p);
mm = wS/S;                    
SPLine = zeros(pwS,pS);        
SPLine(1:mm:end,:) = eye(U);     
for ii = 1:mm-1
    deltexi = ii/mm*dtaui; 
    SPLine(ii+1:mm:end,:) = A + B*deltexi + C*deltexi^2 + D*deltexi^3;
end
wGam0 = SPLine;       
wGam1 = SPLine*xi1;    
end

% dtaui = 2*pi/S;
% taui = (0 : dtaui : 2*p*pi-dtaui)';
% wdtaui = 2*pi/wS;
% wtaui = (0 : wdtaui : 2*p*pi-wdtaui)';
% 
% y = cos(0.5*taui) + sin(3*taui) + cos(taui);
% dy = -0.5*sin(0.5*taui) + 3*cos(3*taui) - sin(taui);
% ddy = -0.5*0.5*cos(0.5*taui) - 9*sin(3*taui) - cos(taui);
% 
% yy = y;
% wyy = SPLine * y;
% dyy = xi1*y;
% wdyy = wGam1*y;
% ddyy = xi2*y;
% wddyy = wGam0*ddyy;
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
% 
% figure(3)
% plot(taui, ddy, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, ddyy, 'r--', 'LineWidth', 1.5);
% plot(wtaui, wddyy, 'b--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','wt');