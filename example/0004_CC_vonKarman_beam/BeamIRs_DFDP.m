function DFDP = BeamIRs_DFDP(~,~,theta,model)
% Derivative with respect to the excitation scale p(1).
DFDP = reshape(model.E*cos(theta(1,:)),size(model.E,1),1,[]);
end
