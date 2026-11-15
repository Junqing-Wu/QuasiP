function [M,C,K,F,E,output] = build_vonKarmanBeamIRs(nElements,forceFormat)
% Reference: "Parallelized computation of quasi-periodic solutions for finite
% element problems: A Fourier series expansion-based shooting method".
% Section 4.2 and Table 3: clamped-pinned von Karman beam.
% M*qdd + C*qd + K*q + F(q) = E*cos(theta), theta = Omega*t.
% [f,J] = F(q) also gives the nonlinear tangent for a single state.
% F(Q) accepts one displacement vector per column.
% Use forceFormat='tensor' for the SSMTool convention F={T2,T3}.
validateattributes(nElements,{'numeric'},{'scalar','real','finite','integer','positive'});
if nargin<2
    forceFormat = 'handle';
end
forceFormat = validatestring(forceFormat,{'handle','tensor'});

%% ====================== Physical parameters ======================
l = 2700;                          % mm
b = 10;
h = 10;
young = 4.5e7;                     % kPa in the kg-mm-s convention
rho = 1780e-9;                     % kg/mm^3
A = b*h;
I = b*h^3/12;
le = l/nElements;
kLinear = 37;                      % N/mm
f1 = 20;                           % N
beta = 2e-4/9;

%% ====================== Assemble the beam ======================
% Element DOFs: [u1,w1,theta1,u2,w2,theta2].
% The quartic stretching energy contains a degree-eight integrand.
% Five Gauss points also integrate the consistent mass exactly.
xi = [-0.906179845938664; -0.538469310105683; 0; ...
       0.538469310105683;  0.906179845938664];
w = (le/2)*[0.236926885056189; 0.478628670499366; ...
            0.568888888888889; 0.478628670499366; 0.236926885056189];
s = (xi+1)/2;
Nu = [1-s,0*s,0*s,s,0*s,0*s];
Bu = repmat([-1,0,0,1,0,0]/le,numel(xi),1);
Nw = transverse_shape(s,le);
Bw = [0*s,(-6*s+6*s.^2)/le,1-4*s+3*s.^2, ...
      0*s,(6*s-6*s.^2)/le,-2*s+3*s.^2];
Bww = [0*s,(-6+12*s)/le^2,(-4+6*s)/le, ...
       0*s,(6-12*s)/le^2,(-2+6*s)/le];
Me = rho*A*(Nu'*(w.*Nu)+Nw'*(w.*Nw));
Ke = young*(A*Bu'*(w.*Bu)+I*Bww'*(w.*Bww));

nFull = 3*(nElements+1);
edof = (1:6)' + 3*(0:nElements-1);
i = repmat(edof,6,1);
j = repelem(edof,6,1);
M = sparse(i(:),j(:),repmat(Me(:),nElements,1),nFull,nFull);
K = sparse(i(:),j(:),repmat(Ke(:),nElements,1),nFull,nFull);

% Clamp [u,w,theta] on the left and pin [u,w] on the right.
free = [4:nFull-3,nFull];
M = M(free,free);
K = K(free,free);
mid = point_shape(0.5,le,edof,nFull);
quarter = point_shape(0.25,le,edof,nFull);
mid = mid(free);
quarter = quarter(free);
K = K+kLinear*(mid*mid');
C = beta*K;
E = f1*mid;
if strcmp(forceFormat,'tensor')
    if isempty(which('sptensor'))
        error('build_vonKarmanBeamIRs:MissingTensorToolbox', ...
            'Add MATLAB Tensor Toolbox to the path for forceFormat=''tensor''.');
    end
    F = nonlinear_tensors(free,edof,young*A*w,Bu,Bw,nFull);
else
    F = @(q) nonlinear_force(q,free,edof,young*A*w,Bu,Bw,nFull);
end
output = [quarter';mid'];          % [w(l/4); w(l/2)] = output*q
end

function Nw = transverse_shape(s,le)
% Cubic Hermite interpolation of transverse displacement.
Nw = [0*s,1-3*s.^2+2*s.^3,le*(s-2*s.^2+s.^3), ...
      0*s,3*s.^2-2*s.^3,le*(-s.^2+s.^3)];
end

function F = nonlinear_tensors(free,edof,EAweights,Bu,Bw,nFull)
% Assemble polynomial coefficients on sorted input indices as in SSMTool.
% The tensor constructor sums repeated monomials and shared element entries.
reduced = zeros(nFull,1);
reduced(free) = 1:numel(free);
F = cell(1,2);
for order = 2:3
    if order==2
        [i,j,k] = ndgrid(1:6);
        subs = [i(:),j(:),k(:)];
        values = EAweights'*(0.5*Bu(:,i(:)).*Bw(:,j(:)).*Bw(:,k(:)) ...
            + Bw(:,i(:)).*Bu(:,j(:)).*Bw(:,k(:)));
    else
        [i,j,k,l] = ndgrid(1:6);
        subs = [i(:),j(:),k(:),l(:)];
        values = 0.5*EAweights'*(Bw(:,i(:)).*Bw(:,j(:)).*Bw(:,k(:)).*Bw(:,l(:)));
    end
    keep = values(:)~=0;
    subs = subs(keep,:);
    subs(:,2:end) = sort(subs(:,2:end),2);
    values = values(keep);
    nTerms = size(subs,1);
    globalSubs = zeros(nTerms*size(edof,2),order+1);
    for ie = 1:size(edof,2)
        idx = edof(:,ie);
        globalSubs((ie-1)*nTerms+(1:nTerms),:) = reduced(idx(subs));
    end
    globalValues = repmat(values(:),size(edof,2),1);
    freeTerms = all(globalSubs>0,2);
    F{order-1} = sptensor(globalSubs(freeTerms,:),globalValues(freeTerms), ...
        repmat(numel(free),1,order+1));
end
end

function shape = point_shape(fraction,le,edof,nFull)
% Evaluate the force or sensor shape at a physical point of the beam.
position = fraction*size(edof,2);
element = floor(position)+1;
Nw = transverse_shape(position-floor(position),le);
shape = sparse(edof(:,element),ones(6,1),Nw(:),nFull,1);
end

function [f,J] = nonlinear_force(q,free,edof,EAweights,Bu,Bw,nFull)
% Assemble the strain operators on the reduced mesh once per discretization.
persistent cache
nElements = size(edof,2);
if isempty(cache) || cache.nElements~=nElements || ~isequal(cache.Bw,Bw) ...
        || ~isequal(cache.EAweights,EAweights)
    ng = numel(EAweights);
    rows = repmat((1:ng)'+ng*(0:nElements-1),6,1);
    cols = repelem(edof,ng,1);
    U = sparse(rows(:),cols(:),repmat(Bu(:),nElements,1),ng*nElements,nFull);
    W = sparse(rows(:),cols(:),repmat(Bw(:),nElements,1),ng*nElements,nFull);
    cache = struct('nElements',nElements,'Bw',Bw,'EAweights',EAweights, ...
        'U',U(:,free),'W',W(:,free),'weights',repmat(EAweights,nElements,1));
end
U = cache.U;
W = cache.W;
weights = cache.weights;
ux = U*q;
wx = W*q;
strainNL = 0.5*wx.^2;
f = U'*(weights.*strainNL)+W'*(weights.*(ux+strainNL).*wx);
if nargout>1
    assert(size(q,2)==1,'Request the tangent for one displacement vector.');
    ng = size(U,1);
    coupling = U'*spdiags(weights.*wx,0,ng,ng)*W;
    J = coupling+coupling.'+W'*spdiags(weights.*(ux+3*strainNL),0,ng,ng)*W;
end
end
