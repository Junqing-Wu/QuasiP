function [xSample,ySample,stabSample,idxSample] = ...
    sample_curve_points(x,y,stability,numSamplePoints,densityWeight)
%SAMPLE_CURVE_POINTS Adaptive sampling along a curve.
%
% densityWeight:
%   0 - Uniform sampling according to curve arc length
%   1 - Uniform sampling according to original point indices
%
% A value between 0.5 and 0.8 is generally recommended.

x = x(:).';
y = y(:).';
stability = stability(:).';

numPoints = numel(x);

if numel(y)~=numPoints || numel(stability)~=numPoints
    error('x, y, and stability must have the same number of elements.');
end

if numSamplePoints<2
    error('numSamplePoints must be at least 2.');
end

if densityWeight<0 || densityWeight>1
    error('densityWeight must be between 0 and 1.');
end

% Normalize the coordinates
xRange = max(x)-min(x);
yRange = max(y)-min(y);

if xRange==0
    xRange = 1;
end

if yRange==0
    yRange = 1;
end

xNormalized = (x-min(x))/xRange;
yNormalized = (y-min(y))/yRange;

% Normalized arc-length increments
arcStep = hypot(diff(xNormalized),diff(yNormalized));

if sum(arcStep)>0
    arcStep = arcStep/sum(arcStep);
else
    arcStep = ones(1,numPoints-1)/(numPoints-1);
end

% Uniform increments according to original point indices
indexStep = ones(1,numPoints-1)/(numPoints-1);

% Combine geometric spacing and original point density
samplingStep = ...
    (1-densityWeight)*arcStep + densityWeight*indexStep;

samplingCoordinate = [0,cumsum(samplingStep)];

% Generate approximately uniformly spaced sampling locations
numSamplePoints = min(numSamplePoints,numPoints);

targetCoordinate = linspace( ...
    samplingCoordinate(1), ...
    samplingCoordinate(end), ...
    numSamplePoints);

idxUniform = interp1( ...
    samplingCoordinate, ...
    1:numPoints, ...
    targetCoordinate, ...
    'nearest');

% Preserve both sides of every stability transition
idxTransition = find(diff(stability)~=0);

idxSample = unique([ ...
    1, ...
    idxUniform, ...
    idxTransition, ...
    idxTransition+1, ...
    numPoints]);

xSample = x(idxSample);
ySample = y(idxSample);
stabSample = stability(idxSample);

end