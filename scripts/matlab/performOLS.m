% performOLS - Function to take a grid, randomly sample some points, and
% atttempt to build a trend using OLS
%
% Syntax:  [h] = performOLS(dat, dim, name)
%
% Inputs:
%    path - matrix containing values
%    dim  - dimensions on grid
%
% Outputs:
%    h - handler to figure
%
% Example:
%    [h] = loadObjects(U, dim, 'Elevation')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Lewis Li
% email: lewisli@stanford.edu
% Website: http://www.
% July 2013; Last revision: 03-July-2013

clear all;
close all;



% Path to true grid
%trueGridPath = ['../data/DS-NonStationary/Reference'];

trueGridPath = ['../data/WLake/WalkerLake'];

% Number of points to use for sampling for reconstructing trend
numHardPoints = 50


% Load grid
[gridDim gridVals gridName] = loadObjects(trueGridPath);

displayRealizations(gridVals, gridDim, gridName);

HDatVals = sampleHardData(gridVals(:,1),gridDim, numHardPoints);
HDatVals(:,1:3) = HDatVals(:,1:3) - 1;

InputVector = [HDatVals(:,1) HDatVals(:,2) ones(numHardPoints,1)];

beta = mvregress(HDatVals(:,1:3), HDatVals(:,4), 'algorithm','ecm');

trendEstimate = zeros(gridDim);
for i = 1:gridDim(1)
   for j = 1:gridDim(2)
       trendEstimate(i,j) = i*beta(1) + j*beta(2) + beta(3);
   end
end

%displayRealizations(trendEstimate, gridDim, 'TrendEstimate');
%