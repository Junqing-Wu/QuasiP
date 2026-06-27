function [Omic,dOmic_dOm] = MT_SHT_Omic(PHI,Om1,Rho,d,Omt,Omnot)
kQ = PHI.kQ;
P  = cell2mat(PHI.P);
lamda = 2*P(Omt)*pi*Rho'*kQ;

v0   = [1,kron(cos(lamda),[1,1])];
v1   = kron(sin(lamda),[0,1]);
Omic = sparse(diag(v0)+diag(-v1,1)+diag(v1,-1));

dOmic_dOm = cell(1,d);
for  ii = 1:d-1
     v0 = [0,zeros(1,2*size(kQ,2))];
     v1 = kron(kQ(ii,:),[0,1]);
     D = diag(v0)+diag(-v1,1)+diag(v1,-1);
     dOmic_dOm{Omnot(ii)} = 2*P(Omt)*pi/Om1*D*Omic;
end
v1 = kron(lamda,[0,1]);
D = diag(v0)+diag(-v1,1)+diag(v1,-1);
dOmic_dOm{Omt} = -1/Om1*D*Omic;
end
