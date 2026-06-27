function dy    = MT_SHT_dy(data,y,theta,Om1,p,N)
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