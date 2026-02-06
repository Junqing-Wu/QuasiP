% 2ND ODEs: M*\ddot{x} + D*\dot{x} + K*x + n(x,\dot{x}) = e  
% 1ST ODEs: B*\dot{X} = A*X + N(X) + E  ----> X = [] 
% B = [D M; M 0]; A = [-K 0; 0 M]; N = [-n(X); 0]; E = [e; 0];

%  \ddot{x} - mu (1-x^2)*\dot{x} + Om0^2*x +alp*x^3 = \sum_{i=1}^d epsd*cos(OMd*t)
%  M =1; + D = -mu; K = Om0^2; n(x,\dot{x}) = mu*x^2*\dot{x} + alp*x^3; e = \sum_{i=1}^d epsd*cos(OMd*t);
%  B = [-mu 1; 1 0]; A = [-Om0^2 0; 0 1]
%  N = [-(mu*X1^2*X2 + alp*X1^3) ; 0]
%  E = [e = \sum_{i=1}^d epsd*cos(OMd*t); 0]

function F_x = VD_F_X(tau,x,p,varargin)

x1 = x(1,:);
x2 = x(2,:);

mu  = p(1);
Om0 = 2;
alp = 0.5;


%%  Ax_x = A;

A = [-Om0^2  0;
       0     1]; 

%%  N_x

N_x = zeros(2,2,numel(x1));

N_x(1,1,:) = - ( 2*mu.*x1.*x2  + 3*alp.* x1.^2 );
N_x(1,2,:) = - ( mu.*x1.^2 );

%%

if nargin == 3 
    F_x = repmat(A,1,1,numel(x1)) + N_x;
else
    F_x.Ax_x   = A;    
    F_x.N_x = N_x;
end
end
