function dy  = MT_SHT_dy(data,y,theta,Om1,p,N)
%MT_SHT_DY  Evaluate the shooting ODE in the scaled phi1 coordinate.
%
% The independent variable is phi1=Om1*t. For first-order models the code
% evaluates dy/dphi1 = f(y,theta,p)/Om1. For second-order models it evaluates
% the equivalent scaled equation using the mass, damping, and stiffness terms.
%
% On input:
%   data : QuasiP model data and ODE flags.
%   y : Current first-order or second-order state.
%   theta : Excitation phase vector.
%   Om1 : Selected integration frequency.
%   p : Model-parameter vector.
%   N : First-order state dimension.
%                          selected frequency, parameters, and state dimension.
%
% On output:
%   dy : Derivative with respect to phi1.
if data.ODE == 1
    x  = y';
else
    x  = y(:,1:N/2)';
    dx = y(:,1+N/2:end)';
end
%------------------------------------------------------
if ~data.isL
    g  = data.f(x,p,theta);
else
    [~,~,A,B,C,~,~,~] = data.l(p);
    if ~data.isABC(1)
        f        = data.f(x,p,theta);
        if ~data.isABC(3)
            C = zeros(N,N);
        end
        g     = B\(f-C*x);
    else
        if ~data.isfwithu
            f        = data.f(x,p,theta);
        else
            f        = data.f(x,dx,p,theta);
        end
        if ~data.isABC(2)
            B = zeros(N/2,N/2);
        end
        if ~data.isABC(3)
            C = zeros(N/2,N/2);
        end
        g     = [dx;A\(f-C*x-B*dx)];
    end
end
dy  =  g'/Om1;
end