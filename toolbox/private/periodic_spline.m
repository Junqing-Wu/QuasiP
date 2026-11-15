function [A,B,C,D] = periodic_spline(U,p) % cubic spline interpolation for periodic solution
%PERIODIC_SPLINE  Construct periodic cubic-spline interpolation coefficients.
%
% Solve the cyclic cubic-spline equations for a periodic grid. The returned
% matrices A, B, C, and D provide the value contribution on each interval:
%   x(delta) = A + B*delta + C*delta^2 + D*delta^3.
% The cyclic endpoint equations make the first and last intervals consistent.
%
% On input:
%   U : Number of periodic grid nodes.
%   p : Period multiplier defining the interval length.
%
% On output:
%   A : Cubic-spline coefficient matrices.
%   B : Cubic-spline coefficient matrices.
%   C : Cubic-spline coefficient matrices.
%   D : Cubic-spline coefficient matrices.
A = eye(U); % Identity matrix (zero-order term, y)
h = 2*p*pi/U; % Uniform grid spacing in [0, 2pπ)

%% Construct linear system for 2nd derivative (periodic boundary)
% C1: Coefficient matrix for 2nd derivative continuity (cyclic tridiagonal)
C1 = diag(ones(U,1)*4*h^2,0) + ...     % Main diagonal: 4h²
    diag(ones(U-1,1)*h^2,1) + ...     % Upper 1st diagonal: h²
    diag(ones(U-1,1)*h^2,-1) + ...    % Lower 1st diagonal: h²
    diag(h^2,U-1) + ...               % Cyclic term: (1,U) position
    diag(h^2,1-U);                    % Cyclic term: (U,1) position

% C2: Right-hand side matrix for 2nd derivative (cyclic boundary)
C2 = diag(-ones(U,1)*6,0) + ...        % Main diagonal: -6
    diag(ones(U-1,1)*3,1) + ...       % Upper 1st diagonal: 3
    diag(ones(U-1,1)*3,-1) + ...      % Lower 1st diagonal: 3
    diag(3,U-1) + ...                 % Cyclic term: (1,U) position
    diag(3,1-U);                      % Cyclic term: (U,1) position

C =  C1\C2; % Solve linear system C1·C = C2 (2nd derivative coefficients)

%% Construct 3rd derivative coefficient matrix (D)
D1 = diag(-ones(U,1)/h/3,0) + ...      % Main diagonal: -1/(3h)
    diag(ones(U-1,1)/h/3,1) + ...     % Upper 1st diagonal: 1/(3h)
    diag(1/h/3,1-U);                  % Cyclic term: (1,U) position
D = D1*C; % 3rd derivative coefficients (from 2nd derivative)

%% Construct 1st derivative coefficient matrix (B)
B1 = diag(-ones(U,1)/h,0) + ...        % Main diagonal: -1/h
    diag(ones(U-1,1)/h,1) + ...       % Upper 1st diagonal: 1/h
    diag(1/h,1-U);                    % Cyclic term: (1,U) position
% Combine terms for 1st derivative (B = B1 - hC - h²D)
B = B1 + (-h)*C + (-h^2)*D;
end