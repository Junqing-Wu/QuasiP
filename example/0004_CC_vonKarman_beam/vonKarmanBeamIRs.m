clearvars;
close all;
clc;

%% ====================== Beam model ======================
% Reference: "Parallelized computation of quasi-periodic solutions for finite
% element problems: A Fourier series expansion-based shooting method".
% Section 4.2.
nElements = 32;
[M,C,K,F,E,output] = build_vonKarmanBeamIRs(nElements);
n = size(M,1);

% The first two undamped frequencies are close to a 1:3 ratio.
omega = sqrt(sort(eig(full(K),full(M))));
fprintf('%d elements, %d DOFs\n',nElements,n);
fprintf('omega_1 = %.6f, omega_2 = %.6f rad/s\n',omega(1),omega(2));

