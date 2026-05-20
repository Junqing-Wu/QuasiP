function PHI = MP_mVCF_COM_PHI(PHI_INITIAL)
PHI = PHI_INITIAL;
d   = PHI.d;
n   = PHI.n;
isCO = false;
%----------------------------------------------------------
U      = cell(1,d);
pS     = cell(1,d);
pwS    = cell(1,d);
%----------------------------------------------------------
xi0    = cell(1,d);
xi1    = cell(1,d);
xi2    = cell(1,d);
%----------------------------------------------------------
wxi1   = cell(1,d);
wxi2   = cell(1,d);
%----------------------------------------------------------
XI0    = 1;
XI1    = cell(1,d);
XI2    = cell(d,d);
%----------------------------------------------------------
Lam1   = cell(1,d);
Lam2   = cell(1,d);
%----------------------------------------------------------
Gam0   = cell(1,d);
Gam1   = cell(1,d);
iGam0  = cell(1,d);
%----------------------------------------------------------
wGam0  = cell(1,d);
wGam1  = cell(1,d);
%----------------------------------------------------------
for ii = 1:d
    ki   = PHI.k{ii};
    pii  = PHI.P{ii};
    Si   = PHI.S{ii};
    wSi  = PHI.S_wave{ii};
    switch PHI.method{ii}
        case{'HB'}
            [xi1{ii},xi2{ii},Gam0{ii},iGam0{ii},...
                wGam0{ii},U{ii},pS{ii},pwS{ii}] = MP_mVCF_ShapeF_HB(ki,pii,Si,wSi);
            xi0{ii}   = eye(U{ii});
            Gam1{ii}  = Gam0{ii}*xi1{ii};
            wGam1{ii} = wGam0{ii}*xi1{ii};
            wxi1{ii}  = xi1{ii};
            wxi2{ii}  = xi2{ii};
        case{'FD'}
            [xi1{ii},xi2{ii},wGam0{ii},...
                wGam1{ii},U{ii},pS{ii},pwS{ii}] = MP_ShapeF_FD(ki,pii,Si,wSi);
            xi0{ii}   = eye(U{ii});
            Gam0{ii}  = eye(U{ii});
            Gam1{ii}  = xi1{ii};
            iGam0{ii} = eye(U{ii});
            wxi1{ii}  = xi1{ii};
            wxi2{ii}  = xi2{ii};
        case{'CO'}
            isCO = true;
            [xi0{ii},xi1{ii},xi2{ii},wxi1{ii},wxi2{ii},wGam0{ii},...
                wGam1{ii},U{ii},pS{ii},pwS{ii}] = MP_ShapeF_CO(ki,pii,Si,wSi);
            Gam0{ii}  = eye(U{ii});
            Gam1{ii}  = wxi1{ii};
            iGam0{ii} = eye(U{ii});
    end
    XI0 = sparse( kron(XI0,xi0{ii}) );
end
%----------------------------------------------------------
PHI.U   = U;
PHI.pS  = pS;
PHI.pwS = pwS;
%----------------------------------------------------------
U = cell2mat(U);
for ii = 1:d
    U1_im1 = prod(U(1:ii-1));  
    Uip1_d = prod(U(ii+1:d));
    nU1_im1 = n*U1_im1;
    Lam1{ii} = sparse( kron( kron(eye(nU1_im1),wxi1{ii}), eye(Uip1_d)) );
    Lam2{ii} = sparse( kron( kron(eye(nU1_im1),wxi2{ii}), eye(Uip1_d)) );
end
%----------------------------------------------------------
for ii =1:d
    XI1{ii} = 1;  
    for jj = 1:d
        if ii == jj
            XI1{ii} = sparse( kron(XI1{ii},xi1{jj}) );
        else
            XI1{ii} = sparse( kron(XI1{ii},xi0{jj}) );            
        end
    end
end
%----------------------------------------------------------
for ii = 1:d
    for jj = 1:d
        XI2{ii,jj} = 1;
        for kk = 1:d
            if kk == ii || kk == jj
                if ii == jj
                   XI2{ii,jj} = sparse( kron(XI2{ii,jj},xi2{kk}) );
                else
                   XI2{ii,jj} = sparse( kron(XI2{ii,jj},xi1{kk}) );
                end
            else
                XI2{ii,jj} = sparse( kron(XI2{ii,jj},xi0{kk}) );
            end
        end
    end
end
%----------------------------------------------------------
S  = cell2mat(PHI.S);
pS = cell2mat(PHI.pS);
p  = cell2mat(PHI.P);
dtau1   = 2*pi/S(1);
Tau1    = 2*p(1)*pi;
tau_1im1 = (0:dtau1:Tau1-dtau1)';
pS_1im1 = pS(1);
for ii = 2:d
    dtaui = 2*pi/S(ii);
    Taui  = 2*p(ii)*pi;
    tau_i = (0:dtaui:Taui-dtaui)';  
    tau_1im1 = [ kron( ones(pS(ii),1),tau_1im1 ) , kron( tau_i,ones(pS_1im1,1) ) ];
    pS_1im1 = pS_1im1*pS(ii);  
end
tau = tau_1im1'; 
%----------------------------------------------------------
wS  = cell2mat(PHI.S_wave);
pwS = cell2mat(PHI.pwS);
p   = cell2mat(PHI.P);
wdtau1   = 2*pi/wS(1);
wTau1    = 2*p(1)*pi;
wtau_1im1 = (0:wdtau1:wTau1-wdtau1)';
pwS_1im1 = pwS(1);
for ii = 2:d
    wdtaui = 2*pi/wS(ii);
    wTaui  = 2*p(ii)*pi;
    wtau_i = (0:wdtaui:wTaui-wdtaui)';  
    wtau_1im1 = [ kron( ones(pwS(ii),1),wtau_1im1 ) , kron( wtau_i,ones(pwS_1im1,1) ) ];
    pwS_1im1 = pwS_1im1*pwS(ii);  
end
wtau = wtau_1im1';
%----------------------------------------------------------
PHI.isCO = isCO;
PHI.tau    = tau;
PHI.wtau   = wtau;
%----------------------------------------------------------
PHI.xi0   = xi0;
PHI.xi1   = xi1;
PHI.xi2   = xi2;
%----------------------------------------------------------
PHI.XI0   = XI0;
PHI.XI1   = XI1;
PHI.XI2   = XI2;
%----------------------------------------------------------
PHI.Lam1  = Lam1;
PHI.Lam2  = Lam2;
%----------------------------------------------------------
PHI.Gam0  = Gam0;
PHI.Gam1  = Gam1;
PHI.iGam0 = iGam0;
%----------------------------------------------------------
PHI.wGam0 = wGam0;
PHI.wGam1 = wGam1;
%----------------------------------------------------------
PHI.Theta = sparse( kron(eye(n),XI0) );
PHI.L     = prod(U);
PHI.nL    = n*PHI.L;
end

