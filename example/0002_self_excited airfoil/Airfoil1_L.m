function [l,lp,A,B,C,Ap,Bp,Cp] = Airfoil1_L(p)
% Parameters of system
miu=12.8;
miubeta=4;
xalpha=0.15;
xbeta=0.2;
ralpha2=0.3;
rbeta2=0.89;
Lbar=0.18;
a=-0.41;
% b=0.118;
ch=0.2;
calpha=0.2;
cbeta=0;
wh=34.6;
walpha=88;
wbeta=60;
% kbeta3=10;
%----------------------------------------------------------
Q = p(1,1);
%----------------------------------------------------------
l  = [1,1,1];
lp = [0,0,1];
%
M11 = miu + miubeta;
M12 = miu*xalpha + miubeta*xbeta - miubeta*Lbar;
M13 = miubeta*xbeta;
M22 = miu*ralpha2 + miubeta*(rbeta2 + Lbar^2) - 2*miubeta*xbeta*Lbar;
M23 = miubeta*rbeta2 - miubeta*xbeta*Lbar;
M33 = miubeta*rbeta2;
A  = [M11 M12 M13; M12 M22 M23; M13 M23 M33];
Ap = [];
%----------------------------------------------------------
B = diag([ch, calpha, cbeta]);
Bp = [];
%----------------------------------------------------------
K11 = miu * (wh/walpha)^2;
K22 = miu*ralpha2 - 2*(Lbar + a)*Q;
K33 = miubeta*rbeta2 * (wbeta/walpha)^2;
K12 =  2*Q;
K21 = 0; K13 = 0; K31 = 0; K23 = 0; K32 = 0;
C = [K11 K12 K13; K21 K22 K23; K31 K32 K33];
%----------------------------------------------------------
dK22_dQ = -2 * (Lbar + a);
dK12_dQ = 2;
Cp = [0, dK12_dQ, 0;
      0, dK22_dQ, 0;
      0, 0      , 0];
%----------------------------------------------------------
end