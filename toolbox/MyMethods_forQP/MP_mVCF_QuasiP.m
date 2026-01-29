function [data, y, J] = MP_mVCF_QuasiP(prob, data, u)

nL = data.nL;
d = data.ode_opts.d;
n = data.DOF;
PHI = data.PHI;
U = cell2mat(PHI.U);
U_prod = prod(U);
S = cell2mat(PHI.S);
S_prod = prod(S);
tau = PHI.tau';

Xd = u(1:nL);        
Xd_n = reshape(Xd,U_prod,n)';
Om = u(nL+1:nL+d);
p  = u(nL+d+1:end);


% 2ND ODEs: M*\ddot{x} + D*\dot{x} + K*x + n(x,\dot{x}) = e  
% 1ST ODEs: B*\dot{X} = A*X + N(X) + E  ----> X = [] 
% B = [D M; M 0]; A = [-K 0; 0 M]; N = [-n(X); 0]; E = [e; 0];
switch data.order
    case{1} % 1st ODEs:

        Xi1 = sparse(U_prod,U_prod);
        for ii = 1:d 
            Xi1 = Xi1 + Om(ii)*PHI.Xi.st{ii};
        end
        Xi0 = PHI.Xi.zth;

        X0 = zeros(n,S_prod);
        parfor nn = 1:n
            X0(nn,:) = U2S_0(PHI,Xd_n(nn,:),U_prod);
        end

        B   = data.b(p);
        B_P = data.bp(p);
        F = data.f(tau,X0,p,'');
        F_X = data.fx(tau,X0,p,'');
        F_P = data.fp(tau,X0,p,'');

        %for y = [ -(B \otimes Xi^1) + A \otimes Xi^0 ]*Xd + [eye \otimes Xi^0]*[Nd + Ed]

        y = - kron(B,Xi1)*Xd;
        if isfield(F, 'A')
            y = y + kron(F.A,Xi0)*Xd;
        end
        if isfield(F, 'N_tau')
            Nd_n = zeros(U_prod,n);
            for nn = 1:n
                is_all_zero = all(F.N_tau(nn,:) == 0);
                if ~is_all_zero
                    Nd_n(:,nn) = S2U_1(PHI,F.N_tau(nn,:),S_prod);
                end
            end
            Nd = reshape(Nd_n,[],1);
            if PHI.isCO1 == 0
               y = y + Nd;
            else 
               y = y + kron(eye(n),Xi0)*Nd;                
            end
        end
        if isfield(F, 'E_tau')
            Ed_n = zeros(U_prod,n);
            for nn = 1:n
                is_all_zero = all(F.E_tau(nn,:) == 0);
                if ~is_all_zero
                    Ed_n(:,nn) = S2U_1(PHI,F.E_tau(nn,:),S_prod);
                end
            end
            Ed = reshape(Ed_n,[],1);
            if PHI.isCO1 == 0
                y = y + Ed;
            else
                y = y + kron(eye(n),Xi0)*Ed;
            end
        end

        % J_Om
        J_Om = zeros(nL,d);
        for ii = 1:d
            J_Om(:,ii) = -kron(B,PHI.Xi.st{ii})*Xd;
        end

        % J_Xd
        J_Xd = - kron(B,Xi1);
        if isfield(F_X, 'Ax_x')
            J_Xd = J_Xd + kron(F_X.Ax_x,Xi0);
        end
        if isfield(F_X, 'N_x')
            Nd_Xd = zeros(nL,nL);
            for ll = 1:n
                for mm = 1:n
                    is_all_zero = all(F_X.N_x(ll,mm,:) == 0);
                    if ~is_all_zero
                        Ndl_Xdm = S2U_2_withX(PHI,F_X.N_x(ll,mm,:),S_prod);
                        Nd_Xd( (ll-1)*U_prod+1:ll*U_prod , ...
                               (mm-1)*U_prod+1:mm*U_prod )  = ...
                        Nd_Xd( (ll-1)*U_prod+1:ll*U_prod , ...
                               (mm-1)*U_prod+1:mm*U_prod ) + Ndl_Xdm;
                    end
                end
            end
            if PHI.isCO1 == 0
                J_Xd = J_Xd + Nd_Xd;
            else
                J_Xd = J_Xd + kron(eye(n),Xi0)*Nd_Xd;
            end
        end

        % J_p
        J_p = zeros(nL,numel(p));

        for ii = 1:numel(p)
            B_Pii = B_P(:,:,ii);
            if any(B_Pii(:) ~= 0) 
                J_p(:,ii) = J_p(:,ii) + kron(-B_Pii,Xi1)*Xd;
            end
            if isfield(F_P, 'A_p')
                A_Pii = F_P.A_p(:,:,ii);
                if any(A_Pii(:) ~= 0)
                    J_p(:,ii) = J_p(:,ii) + kron(A_Pii,Xi0)*Xd;
                end
            end
            if isfield(F_P,'NE_p')
                NE0_Pii = reshape(F_P.NE_p(:,ii,:),n,[]);
                if any(NE0_Pii(:) ~= 0)
                    NEd_Pii = zeros(U_prod,n);
                    for nn = 1:n
                        is_all_zero = all(NE0_Pii(nn,:) == 0);
                        if ~is_all_zero
                            NEd_Pii(:,nn) = S2U_1(PHI,NE0_Pii(nn,:),S_prod);
                        end
                    end
                    NEd_Pii = reshape(NEd_Pii,[],1);
                    if PHI.isCO1 == 0
                        J_p(:,ii) = J_p(:,ii) + NEd_Pii;
                    else
                        J_p(:,ii) = J_p(:,ii) + kron(eye(n),Xi0)*NEd_Pii;
                    end
                end
            end

        end

        J = [J_Xd,J_Om,J_p];
    case{2}
end

end

function X0 = U2S_0(PHI,Xd,U_prod)

Siip1Td = 1;
U1Tii = U_prod;
U = cell2mat(PHI.U);
S = cell2mat(PHI.S);

Xii = Xd;
for kk = PHI.d:-1:1
    U1Tii = U1Tii/U(kk);
    Xii_wave = pagetranspose( reshape(Xii,Siip1Td,U(kk),U1Tii) );
    Xii_wave = reshape(Xii_wave,U(kk),Siip1Td*U1Tii);
    switch PHI.method{kk}
        case{'HB'}
            GAM = PHI.Gamma_i{1,kk};
            Xiiminus1_bar = GAM*Xii_wave;
        case{'CO'}
            Xiiminus1_bar = Xii_wave;
        case{'FD'}
            Xiiminus1_bar = Xii_wave;
    end
    Siip1Td = Siip1Td*S(kk);
    Xii = reshape( Xiiminus1_bar, Siip1Td, U1Tii  );
end
X0 = Xii';

end


function Nd = S2U_1(PHI,N0,S_prod)

U = cell2mat(PHI.U);
S = cell2mat(PHI.S);

S_iiTd = S_prod;
U_1Tiiminus1 = 1;

Niiminus1 = N0;

for kk = 1:PHI.d
    S_iiTd = S_iiTd/S(kk);
    Niiminus1_bar = pagetranspose(reshape(Niiminus1,U_1Tiiminus1,S(kk),S_iiTd) );
    Niiminus1_bar = reshape(Niiminus1_bar,S(kk),U_1Tiiminus1*S_iiTd);

    switch PHI.method{kk}
        case{'HB'}
            iGAM = PHI.Gamma_i{3,kk};            
            Fnlii = iGAM*Niiminus1_bar;
        case{'CO'}
            Fnlii = Niiminus1_bar;
        case{'FD'}
            Fnlii = Niiminus1_bar;
    end
    U_1Tiiminus1 = U_1Tiiminus1*U(kk);
    Niiminus1 = reshape(Fnlii,U_1Tiiminus1,S_iiTd);
end
Nd = Niiminus1;
end

function Nd_Xd = S2U_2_withX(PHI,N0_X0,S_prod)
U = cell2mat(PHI.U);
S = cell2mat(PHI.S);
S_iiTd = S_prod;
U_1Tiiminus1 = 1;

N_Xiiminus1 = N0_X0;

for kk = 1:PHI.d
    S_iiTd = S_iiTd/S(kk);

    N_Xiiminus1_bar = pagetranspose( reshape(N_Xiiminus1,U_1Tiiminus1^2,S(kk),S_iiTd) );
    N_Xiiminus1_bar = reshape( N_Xiiminus1_bar, S(kk),1,U_1Tiiminus1^2*S_iiTd );
    
    switch PHI.method{kk}
        case{'HB'}
            GAM = PHI.Gamma_i{1,kk};
            iGAM = PHI.Gamma_i{3,kk}; 
            dFnl_dZiiminus1_bar_SU = repmat(N_Xiiminus1_bar,1,U(kk));
            dFnl_dZii = dFnl_dZiiminus1_bar_SU.*GAM;
            dFnl_dZii = pagemtimes( iGAM,dFnl_dZii );
        case{'CO'}
            dFnl_dZiiminus1_bar_SS = zeros(S(kk),S(kk),U_1Tiiminus1^2*S_iiTd);
            for jj = 1:U_1Tiiminus1^2*S_iiTd
                dFnl_dZiiminus1_bar_SS(:,:,jj) = diag(N_Xiiminus1_bar(:,:,jj));
            end
            dFnl_dZii = dFnl_dZiiminus1_bar_SS;
        case{'FD'}
            dFnl_dZiiminus1_bar_SS = zeros(S(kk),S(kk),U_1Tiiminus1^2*S_iiTd);
            for jj = 1:U_1Tiiminus1^2*S_iiTd
                dFnl_dZiiminus1_bar_SS(:,:,jj) = diag(N_Xiiminus1_bar(:,:,jj));
            end
            dFnl_dZii = dFnl_dZiiminus1_bar_SS;
    end
    N_Xiiminus1 = pagetranspose( reshape(dFnl_dZii,U(kk),U(kk)*U_1Tiiminus1,U_1Tiiminus1,S_iiTd) );
    U_1Tiiminus1 = U_1Tiiminus1*U(kk);
    N_Xiiminus1 = reshape(N_Xiiminus1,U_1Tiiminus1^2,S_iiTd);
end
Nd_Xd = transpose( reshape(N_Xiiminus1,U_1Tiiminus1,U_1Tiiminus1) );
end