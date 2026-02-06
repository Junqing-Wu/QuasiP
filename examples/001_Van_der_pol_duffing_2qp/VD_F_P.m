% 2ND ODEs: M*\ddot{x} + D*\dot{x} + K*x + n(x,\dot{x}) = e  
% 1ST ODEs: B*\dot{X} = A*X + N(X) + E  ----> X = [] 
% B = [D M; M 0]; A = [-K 0; 0 M]; N = [-n(X); 0]; E = [e; 0];

%  \ddot{x} - mu (1-x^2)*\dot{x} + Om0^2*x +alp*x^3 = \sum_{i=1}^d epsd*cos(OMd*t)
%  M =1; + D = -mu; K = Om0^2; n(x,\dot{x}) = mu*x^2*\dot{x} + alp*x^3; e = \sum_{i=1}^d epsd*cos(OMd*t);
%  B = [-mu 1; 1 0]; A = [-Om0^2 0; 0 1]
%  N = [-(mu*X1^2*X2 + alp*X1^3) ; 0]
%  E = [e = \sum_{i=1}^d epsd*cos(OMd*t); 0]

function F_p = VD_F_P(tau,x,p,varargin)

% F(tau_1^e,x,p) = A(p)x + N(x,p) + E(tau_1^e,p)
% A(p)x         ---> linear part
% N(x,p)        ---> nonlinear part
% E(tau_1^e,p)  ---> excitation part

tau1 = tau(1,:);
tau2 = tau(2,:); 

x1 = x(1,:);
x2 = x(2,:);

mu  = p(1);
Om0 = p(2);
alp = p(3);
eps1 = p(4);
eps2 = p(5);

%% A_p
% A = [-K 0; 0 M]

A_p = zeros(2,2,numel(p));

A_p(:,:,2) = [-2*Om0  0;
                0     0]; 

%% N_p

N_p = zeros(2,numel(p),numel(x1));

N_p(1,1,:) = - ( x1.^2.*x2 );
N_p(1,3,:) = - ( x1.^3 );

%% E_p

E_p = zeros(2,numel(p),numel(x1));

E_p(1,4,:) = cos(tau1) ;
E_p(1,5,:) = cos(tau2) ;

%% 

if nargin == 3 
    A_p = pagemtimes(A_p,x);
    A_p = permute(A_p,[1,3,2]);

    F_p = A_p + N_p + E_p;
else
    F_p.A_p = A_p;
    F_p.NE_p = N_p + E_p;
end

end