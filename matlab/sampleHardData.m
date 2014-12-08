% sampleHardData - Samples a set of hard data points from a exhaustive grid
%
% Syntax:  [hDatVals] = sampleHardData(grid, gridDim, numPts)
%
% Inputs:
%     grid  - exhaustive grid
%  gridDim  - exhaustive grid dimensions
% numPoints - number of points we want to sample
%
% Outputs:
%    hDatVals - numPtsx4 matrix containing sampled hard data coordinates
%    and values
%
% Example:
%    [hDatVals] = sampleHardData(DEM, DEMDim, 200)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Lewis Li
% email: lewisli@stanford.edu
% Website: Website: http://scrf.stanford.edu/
% July 2013; Last revision: 12-July-2013

function[HDatVals] = sampleHardData(grid, gridDim, numPoints)

% Set random seed [optional]
rng(18);

% Generate random indices
randIndex = rand(numPoints,2);
randIndex(:,1) = ceil(randIndex(:,1)*gridDim(1));
randIndex(:,2) = ceil(randIndex(:,2)*gridDim(2));
randIndex = [randIndex ones(numPoints,2)];

% Reshape true grid
trueU = reshape(grid,gridDim);

% Write points
for i = 1:numPoints
    randIndex(i,4) = trueU(randIndex(i,1),randIndex(i,2));
end

HDatVals = randIndex;

end