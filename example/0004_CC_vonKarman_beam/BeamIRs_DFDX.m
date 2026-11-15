function DFDX = BeamIRs_DFDX(x,~,~,model)
% QuasiP expects one tangent matrix for each sampled state.
if size(x,2)==1
    [~,J] = model.fnl(x);
    DFDX = -J;
    return
end
n = size(x,1);
DFDX = zeros(n,n,size(x,2));
for ii = 1:size(x,2)
    [~,J] = model.fnl(x(:,ii));
    DFDX(:,:,ii) = -J;
end
end

