function [l,lp,A,B,C,Ap,Bp,Cp] = BeamIRs_L(~,model)
% Constant mass, damping and stiffness in the QuasiP convention.
l = [1,1,1];
lp = [0,0,0];
A = model.M;
B = model.C;
C = model.K;
Ap = [];
Bp = [];
Cp = [];
end


