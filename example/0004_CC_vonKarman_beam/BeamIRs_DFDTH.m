function DFDTH = BeamIRs_DFDTH(~,p,theta,model)
% Derivative with respect to the prescribed excitation phase.
DFDTH = reshape(-model.E*(p(1,:).*sin(theta(1,:))),size(model.E,1),1,[]);
end
