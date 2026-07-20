function [data, y, J] = MT_SHT_HT_QuasiP_PC(prob, data, u)
PHI   = data.PHI;
NL    = PHI.nL;
N     = PHI.N;
d     = data.ode_opts.d;
e     = data.ode_opts.e;
a     = data.pdim;
%----------------------------------------------------------
y1   = zeros(d-e,1);        % 未知基频的相位条件PC的残差
Jx1  = zeros(d-e,NL);       % PC关于Xd的导数
Jom1 = zeros(d-e,d);        % PC关于Om的导数
Jp1  = zeros(d-e,a);        % PC关于p 的导数
%----------------------------------------------------------
ys  = reshape(u(1:NL),[],N);
Om  = u(1+NL:d+NL);
p = u(d+NL+1:end);
Om1 = Om(PHI.Omt);
%============================================================
%>>>>>>>> (I_N \otimes iG^\top*iG)dys1   垂直 ys1r <<<<<<<<<<
%>>>>>>>> (I_N \otimes iG^\top*Xi*iG)ys1 垂直 ys1r <<<<<<<<<<
%============================================================
if d > e
    kk = e + 1;
    loca0   = PHI.local_HT{1};
    ys1     = ys(loca0,:);
    Loca0 = PHI.Local_HT{1};
    for ii = 1:d-e
        if kk == PHI.Omt
            theta1  = PHI.Varphi(loca0,:)';
            if PHI.MATLABsolver
                dys1    = MT_SHT_dy(data,ys1,theta1,Om1,repmat(p,1,PHI.pS2d),N);
            else
                ys1(:,N/2+1:end)  = ys1(:,N/2+1:end)*Om1;
                dys1    = MT_SHT_dy(data,ys1,theta1,Om1,repmat(p,1,PHI.pS2d),N);
                dys1(:,N/2+1:end) = dys1(:,N/2+1:end)/Om1;
            end
            Jx1(ii,Loca0) = kron(eye(N),PHI.iGTiG)*reshape(dys1,[],1);
        else 
            Jx1(ii,Loca0) = kron(eye(N),PHI.iGTXI1iG{kk})*reshape(ys1,[],1);
        end
        kk = kk+1;
    end
end
%==========================================================
y = y1;
J = [Jx1, Jom1, Jp1];
end
