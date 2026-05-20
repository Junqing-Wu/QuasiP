function f = my_dxdt(t,x,p,OM,F,L,isfwithu)
f = zeros(size(x));
if isempty(OM)
    theta = [];
else
    theta = OM*t;
end
if isempty(L)
    f = F(x,p,theta);
else
    [l,~,A,B,C] = L(p);
    if l(1) == 1
        n = length(x)/2;
        f(1:n) = x(n+1:end,:);
        if l(2) == 1
            f(1+n:end) = f(1+n:end) - B*x(n+1:end,:);
        end
        if l(3) == 1
            f(1+n:end) = f(1+n:end) - C*x(1:n);
        end
        if isfwithu
            ne = F(x(1:n),x(n+1:end,:),p,theta);
        else
            ne = F(x(1:n),p,theta);
        end
        f(1+n:end) = f(1+n:end) + ne;
        f(1+n:end) = A\f(1+n:end);
    else
        if l(3) == 1
            f = f - C*x;
        end
        ne = F(x,p,theta);
        f = f + ne;
        f = B\f;
    end
end
end