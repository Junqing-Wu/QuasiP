function Cons = MT_SHT_get_MDK_Cons(data, p, n, a, Om1)
[~,~,A,B,C,Ap,Bp,Cp] = data.l(p);
%-------------------------------------------------------------------------------------
if data.isABC(2) ~= 1
    B = zeros(n,n);end
if data.isABC(3) ~= 1
    C = zeros(n,n);end
if data.isABCwithp(1) ~= 1
    Ap = zeros(n,n,a);end
if data.isABCwithp(2) ~= 1
    Bp = zeros(n,n,a);end
if data.isABCwithp(3) ~= 1
    Cp = zeros(n,n,a);end
%-------------------------------------------------------------------------------------
Cons.A         = A;
Cons.Aw        = A*Om1^2;
Cons.invA      = inv(A);
%-------------------------------------------------------------------------------------
Cons.B         = B;
Cons.Bw        = B*Om1;
%-------------------------------------------------------------------------------------
Cons.C         = C;
%-------------------------------------------------------------------------------------
if data.isp
    Cons.A_p     = Ap;
    Cons.B_p     = Bp;
    Cons.C_p     = Cp;
else
    Cons.A_p     = zeros(n,n,a);
    Cons.B_p     = zeros(n,n,a);
    Cons.C_p     = zeros(n,n,a);
end
end
