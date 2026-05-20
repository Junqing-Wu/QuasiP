clc
clear

% Verification Program for Cubic Spline Interpolation
S = 2^5;    % Coarse grid intervals per 2π
wS = 2^8;   % Fine grid intervals for interpolation
p  = 2;     % Period scaling parameter (integer)

pS = p*S;       % Total number of coarse intervals
pwS = p*wS;     % Total number of fine intervals

n_node = pS + 1;
mm = wS/S;
n_interp = pwS + 1;

% Generate grid points
dtau = 2*pi/S;
taui = (0 : dtau : 2*p*pi)';
wdtau = 2*pi/wS;
wtaui = (0 : wdtau : 2*p*pi)';

% Test function (non-periodic)
y = cos(0.5*taui) + sin(3*taui) + cos(sin(taui)) + 0.1*taui;

figure(1)
plot(taui, y, 'g-', 'LineWidth', 1.5);
hold on;

[A,B,C,D] = cubic_spline(pS, p);

SPLine = zeros(n_interp, n_node);
SPLine(1:mm:end, :) = eye(n_node);

for seg = 1:pS
    for ii = 1:mm-1
        deltexi = ii/mm * dtau;
        idx = (seg-1)*mm + ii + 1;
        SPLine(idx, :) = A(seg, :) + B(seg, :)*deltexi + C(seg, :)*deltexi^2 + D(seg, :)*deltexi^3;
    end
end

Y = SPLine * y;
plot(wtaui, Y, 'r--', 'LineWidth', 1.5);
hold off;
grid on;