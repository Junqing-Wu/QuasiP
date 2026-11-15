function F = BeamIRs_F(x,p,theta,model)
% p(1) scales the nominal 20 N excitation.
F = model.E*(p(1,:).*cos(theta(1,:))) - model.fnl(x);
end
