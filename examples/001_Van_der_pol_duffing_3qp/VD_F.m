% 2ND ODEs: M*\ddot{x} + D*\dot{x} + K*x + n(x,\dot{x}) = e  
% 1ST ODEs: B*\dot{X} = A*X + N(X) + E  ----> X = [] 
% B = [D M; M 0]; A = [-K 0; 0 M]; N = [-n(X); 0]; E = [e; 0];

%  \ddot{x} - mu (1-x^2)*\dot{x} + Om0^2*x +alp*x^3 = \sum_{i=1}^d epsd*cos(OMd*t)
%  M =1; D = -mu; K = Om0^2; n(x,\dot{x}) = mu*x^2*\dot{x} + alp*x^3; e = \sum_{i=1}^d epsd*cos(OMd*t);
%  B = [-mu 1; 1 0]; 
%  A = [-Om0^2 0; 0 1]
%  N = [-(mu*X1^2*X2 + alp*X1^3) ; 0]
%  E = [e = \sum_{i=1}^d epsd*cos(OMd*t); 0]


function F = VD_F(tau,x,p,varargin)

% F(tau_1^e,x,p) = A(p)x + N(x,p) + E(tau_1^e,p)
% A(p)x         ---> linear part
% N(x,p)        ---> nonlinear part
% E(tau_1^e,p)  ---> excitation part

tau1 = tau(1,:);
tau2 = tau(2,:); 
tau3 = tau(3,:); 

x1 = x(1,:);
x2 = x(2,:);

mu  = p(1);
Om0 = 2;
alp = 0.5;
eps1 = 2;
eps2 = 1;

%% A(p) 
% A = [-K 0; 0 M]

A = [-Om0^2  0;
       0     1]; 

%% N(x,p)

N_tau = zeros(2,numel(x1));
N_tau(1,:) = - ( mu.*x1.^2.*x2  + alp.* x1.^3 );

%% E(tau_1^e,p)

E_tau = zeros(2,numel(x1));
E_tau(1,:) = eps1.*cos(tau1) + eps2.*cos(tau2) + 0.5.*cos(tau3) ;

E_MF.harmonic_order = [1,0,0;     % 1*tau1+0*tau2
                       0,1,0      % 0*tau1+1*tau2
                       0,0,1];    
E_MF.amplitude      = [eps1,0;  % [eps1;0]*cos(tau1)   % a*sin(tau1)  eps1-1i*a
                       eps2,0;  % [eps2;0]*cos(tau2) 
                       0.5,0];  % [0.5;0]*cos(tau3)
                        

%% 

if nargin == 3 
    F = A.*x + N_tau + E_tau;
elseif nargin == 4
    F.A = A;
    F.N_tau = N_tau;
    F.E_tau = E_tau; 
else
    F.A = A;
    F.E_MF = E_MF;     
end

end