function [xi0,xi1,xi2,wxi1,wxi2,wGam0,wGam1,U,pS,pwS] = MP_ShapeF_CO(m,p,S,wS)
tol = 1e-6;
if abs(mod(wS, S)) > tol
    warning('S_wave is not an integer multiple of S!');
end
pS  = p*S;
U   = p*S;
pwS = p*wS;
N = U/m;
tol = 1e-6;
if abs(mod(U, m)) > tol
    warning('number of sub-interval in CO is not an integer!');
end
%----------------------------------------------------------
dtaui = 2*pi/S;
%----------------------------------------------------------
c_node = gauss_node(m);
c_node = m/2*(c_node+1);
base_p = (0:m-1)';
%----------------------------------------------------------
a1 = repmat( (0:m)',1,m+1 );             
a1 = a1(:);                              
n = m+1;
a1(1:n+1:n*n) = [];                      
a1 = reshape(a1,n-1,n);                  
a2 = repmat( (0:m),m,1 );                
%----------------------------------------------------------
lo_denominator = prod( a2-a1 );          
molecule_w = kron(c_node',ones(m,m+1)) - kron( ones(1,m),a1 );  
molecule_b = kron(base_p',ones(m,m+1)) - kron( ones(1,m),a1 );  
%----------------------------------------------------------
l0_molecule_w = prod(molecule_w);                            
l0_w = reshape( l0_molecule_w./repmat( lo_denominator,1,m ), m+1, m)';  
%----------------------------------------------------------
l1_molecule_w = zeros(size(l0_molecule_w)); 
l1_molecule_b = zeros(size(l0_molecule_w));
for ii = 1:m
    a4 = molecule_w;
    a4(ii,:) = [];
    l1_molecule_w = prod(a4) + l1_molecule_w;  

    a5 = molecule_b;
    a5(ii,:) = [];
    l1_molecule_b = prod(a5) + l1_molecule_b;  
end
%----------------------------------------------------------
l1_w = reshape( l1_molecule_w./repmat( lo_denominator,1,m ), m+1, m)'/dtaui; 
l1_b = reshape( l1_molecule_b./repmat( lo_denominator,1,m ), m+1, m)'/dtaui; 
%----------------------------------------------------------
Lag0_w_a = kron(eye(N),l0_w(:,1:m));     
Lag0_w_b = kron(eye(N),[zeros(m,m-1),l0_w(:,1+m)]); 
Lag0_w_b = [Lag0_w_b(:,end),Lag0_w_b(:,1:end-1)];     
Lag0_w = Lag0_w_a + Lag0_w_b;
%----------------------------------------------------------
Lag1_w_a = kron(eye(N),l1_w(:,1:m));     
Lag1_w_b = kron(eye(N),[zeros(m,m-1),l1_w(:,1+m)]);  
Lag1_w_b = [Lag1_w_b(:,end),Lag1_w_b(:,1:end-1)];     
Lag1_w = Lag1_w_a + Lag1_w_b;            
%----------------------------------------------------------
Lag1_b_a = kron(eye(N),l1_b(:,1:m));     
Lag1_b_b = kron(eye(N),[zeros(m,m-1),l1_b(:,1+m)]);  
Lag1_b_b = [Lag1_b_b(:,end),Lag1_b_b(:,1:end-1)];     
Lag1_b = Lag1_b_a + Lag1_b_b;
%----------------------------------------------------------
xi0 = Lag0_w;
xi1 = Lag1_w; 
xi2 = Lag1_w*Lag1_b;
wxi1 = Lag1_b;        
wxi2 = Lag1_b*Lag1_b;
%----------------------------------------------------------
mm_hr = pwS/N;                     
c_node_hr = (0:mm_hr-1)'/mm_hr*m; 
%----------------------------------------------------------
a1 = repmat( (0:m)',1,m+1 );
a1 = a1(:);
n = m+1;
a1(1:n+1:n*n) = [];
a1 = reshape(a1,n-1,n);
a2 = repmat( (0:m),m,1 );
lo_denominator = prod( a2-a1 );
molecule_w_hr = kron(c_node_hr',ones(m,m+1)) - kron( ones(1,mm_hr),a1 );  
l0_molecule_w_hr = prod(molecule_w_hr);                                  
l0_w_hr = reshape( l0_molecule_w_hr./repmat( lo_denominator,1,mm_hr ), m+1, mm_hr)';  
l1_molecule_w_hr = zeros(1,size(molecule_w_hr,2));
for ii = 1:m
    a4 = molecule_w_hr;
    a4(ii,:) = [];
    l1_molecule_w_hr = prod(a4) + l1_molecule_w_hr;  
end
l1_w_hr = reshape( l1_molecule_w_hr./repmat( lo_denominator,1,mm_hr ), m+1, mm_hr)'/dtaui;
%----------------------------------------------------------
Lag0_w_a_hr = kron(eye(N),l0_w_hr(:,1:m));                     
Lag0_w_b_hr = kron(eye(N),[zeros(mm_hr,m-1),l0_w_hr(:,1+m)]);  
Lag0_w_b_hr = [Lag0_w_b_hr(:,end),Lag0_w_b_hr(:,1:end-1)];     
wGam0 = Lag0_w_a_hr + Lag0_w_b_hr; 
%----------------------------------------------------------
Lag1_w_a_hr = kron(eye(N),l1_w_hr(:,1:m));                     
Lag1_w_b_hr = kron(eye(N),[zeros(mm_hr,m-1),l1_w_hr(:,1+m)]);  
Lag1_w_b_hr = [Lag1_w_b_hr(:,end),Lag1_w_b_hr(:,1:end-1)];     
wGam1 = Lag1_w_a_hr + Lag1_w_b_hr;                         
end

%% Subfunction: Generate Gauss collocation nodes (symbolic calculation)
function c_node = gauss_node(m)
% GAUSS_NODE - Generate Gauss collocation nodes via symbolic differentiation
%
% INPUT:
%   m    - Degree of the collocation polynomial
%
% OUTPUT:
%   c_node - Sorted Gauss collocation nodes in [-1, 1] (roots of m-th Legendre polynomial)
%
% DESCRIPTION:
%   Calculates Gauss nodes by solving the m-th derivative of the Legendre polynomial ((t^2-1)^m),
%   which are the optimal collocation points for high-precision numerical integration/derivation.
g_l = ((sym('t'))^2-1)^m;       % Legendre polynomial (symbolic expression)
dg_l = diff(g_l,m);             % m-th derivative (Gauss node condition: roots of m-th Legendre polynomial)
c_node = double(solve(dg_l));   % Solve for symbolic roots -> numerical Gauss nodes
c_node = sort(c_node);          % Sort nodes in ascending order [-1, 1]
end


% dtaui = 2*pi/S;
% taui = (0 : dtaui : 2*p*pi-dtaui)';
% wdtaui = 2*pi/wS;
% wtaui = (0 : wdtaui : 2*p*pi-wdtaui)';
% 
% rDtaui = 2*pi/(S/m);
% rPtau  = dtaui*c_node;
% rtaui = (0 : rDtaui : 2*p*pi-rDtaui)';
% rtaui = kron(rtaui,ones(m,1)) + kron(ones(N,1),rPtau);
% 
% y = cos(0.5*taui) + sin(3*taui) + cos(taui);
% dy = -0.5*sin(0.5*taui) + 3*cos(3*taui) - sin(taui);
% ddy = -0.5*0.5*cos(0.5*taui) - 9*sin(3*taui) - cos(taui);
% 
% yy = y;
% ryy = xi0 * y;
% wyy = wGam0 * y;
% 
% dyy = wxi1*y;
% rdyy = xi1*y;
% wdyy = wGam1 * y;
% 
% ddyy = wxi2*y;
% rddyy = xi2*y;
% wddyy = wGam1 * dy;
% 
% figure(1)
% plot(taui, y, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, yy, 'r--', 'LineWidth', 1.5);
% plot(rtaui, ryy, 'b--', 'LineWidth', 1.5);
% plot(wtaui, wyy, 'k--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','rt','wt');
% 
% figure(2)
% plot(taui, dy, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, dyy, 'r--', 'LineWidth', 1.5);
% plot(rtaui, rdyy, 'b--', 'LineWidth', 1.5);
% plot(wtaui, wdyy, 'k--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','rt','wt');
% 
% figure(3)
% plot(taui, ddy, 'g-', 'LineWidth', 1.5);
% hold on;
% plot(taui, ddyy, 'r--', 'LineWidth', 1.5);
% plot(rtaui, rddyy, 'b--', 'LineWidth', 1.5);
% plot(wtaui, wddyy, 'k--', 'LineWidth', 1.5);
% hold off;
% grid on;
% legend('t','tt','rt','wt');
