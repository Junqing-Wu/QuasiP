function [M,D,K,F] = MT_SHT_MDKF_var(data,y,phi1,phi,Om1,Om2d,Omt,Omnot,Rho,p,Cons,n,d,e,a)
% 
q = y(1:n);
u = y(n+1:end);
v = u*Om1;
theta = phi;
theta(Omt) = phi1;
theta(Omnot) = phi(Omnot) + Rho*phi1;
%-------------------------------------------------------------
if data.isfwithu
    f           = data.f(q,v,p,theta);
    [dfdq,dfdv] = data.fx(q,v,p,theta);
    if isempty(dfdv)
        dfdv = zeros(n,n);
    end
    if e ~= 0
        dfdth = zeros(n,d);
        dfdth(:,1:e) = data.fth(q,v,p,theta);
    else
        dfdth = zeros(n,d);
    end
    if data.isfwithp
        dfdp  = data.fp(q,v,p,theta);
    else
        dfdp  = zeros(n,a);
    end
else
    f    = data.f(q,p,theta);    
    dfdq = data.fx(q,p,theta);
    dfdv = zeros(n,n);
    if e ~= 0
        dfdth = zeros(n,d);
        dfdth(:,1:e) = data.fth(q,p,theta);
    else
        dfdth = zeros(n,d);
    end
    if data.isfwithp
        dfdp  = data.fp(q,p,theta);
    else
        dfdp  = zeros(n,a);
    end
end
%-------------------------------------------------------------
M = (Cons.Aw);
D = (Cons.Bw-Om1*dfdv);
K = (Cons.C-dfdq);
%-------------------------------------------------------------
F_qqs = zeros(n,n);
F_qus = zeros(n,n);
F_qOm = zeros(n,d);

F_qOm(:,Omt) =  - 2/Om1 *(f - Cons.C*q - Cons.Bw*u) - Cons.B*u ; 
for ii = 1:d-1
    F_qOm(:,Omt)       = F_qOm(:,Omt) - phi1*Om2d(ii)/Om1^2*dfdth(:,Omnot(ii));
    F_qOm(:,Omnot(ii)) = phi1/Om1*dfdth(:,Omnot(ii));
end
F_qp  = zeros(n,a);
for aa = 1:a
    F_qp(:,aa) = F_qp(:,aa) - Cons.A_p(:,:,aa)*Cons.invA* (f - Cons.C*q + Cons.Bw*u) +  ...
                              ( dfdp(:,aa) - Cons.C_p(:,:,aa)*q - Om1*Cons.B_p(:,:,aa)*u ); 
end
F = [F_qqs,F_qus,F_qOm,F_qp];
%-------------------------------------------------------------
end