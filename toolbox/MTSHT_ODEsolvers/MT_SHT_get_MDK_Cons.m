function Cons = MT_SHT_get_MDK_Cons(data, p, Om1, n, a)
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
Cons.Aw        = Om1^2*A;
invAw          = inv(Cons.Aw);
Cons.invAw_Om1 = -2/Om1^3.*invAw;
%-------------------------------------------------------------------------------------
Cons.Bw        = Om1*B;
Cons.Bw_Om1    = B;
%-------------------------------------------------------------------------------------
Cons.C         = C;
%-------------------------------------------------------------------------------------
if data.isp
    Cons.invAw_p = pagemtimes( pagemtimes( -invAw, Om1^2*Ap ), invAw );

    Cons.Bw_p    = Om1*Bp;
    Cons.C_p     = Cp;
else
    Cons.invAw_p = zeros(n,n,a);
    Cons.Bw_p    = zeros(n,n,a);
    Cons.C_p     = zeros(n,n,a);
end
end
