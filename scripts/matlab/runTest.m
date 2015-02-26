% samplePoints.m
%
% Sample points from Grid
%
% Author: Lewis Li
% email: lewisli@stanford.edu
% Website: Website: http://scrf.stanford.edu/
% July 2013; Last revision: 12-July-2013

% Clear workspace
close all;
clear all;
% ------------- ESTIMATION PARAMETERS --------------
trueDEMPath = '../data/WLakeRaw/';
trueDEM = 'WalkerLakeRaw';
% --------------------------------------------------
% --------------------------------------------------
% [Option 1: We randomly sample the true grid for hard data]
% --------------------------------------------------
[TruthDim TruthVals TruthNames] = loadObjects([trueDEMPath trueDEM]);

numPts = [25 50 100 200];

for i = 1:length(numPts)
    
    numHardPoints = numPts(i);
    HDatVals = sampleHardData(TruthVals(:,1),TruthDim, numHardPoints);
    HDatVals(:,1:3) = HDatVals(:,1:3) - 1;
    
    outputFileName = [trueDEMPath 'hardData' num2str(numHardPoints)];
    
    
    fid=fopen(outputFileName, 'w');
    txt = [ 'random samples\n4\nX\nY\nZ\nU-random-sampling\n'];
    S = sprintf(txt);
    fprintf(fid,'%s',S);
    

    dlmwrite([outputFileName], HDatVals, 'delimiter', ' ', 'precision', 6,'-append');
    
end
numPts = [25 50 100 200];

for i = 1:length(numPts)
    
    numHardPoints = numPts(i);
    HDatVals = sampleHardData(TruthVals(:,1),TruthDim, numHardPoints);
    HDatVals(:,1:3) = HDatVals(:,1:3) - 1;
    
    outputFileName = [trueDEMPath 'hardData' num2str(numHardPoints)];
    
    
    fid=fopen(outputFileName, 'w');
    txt = [ 'random samples\n4\nX\nY\nZ\nU-random-sampling\n'];
    S = sprintf(txt);
    fprintf(fid,'%s',S);
    

    dlmwrite([outputFileName], HDatVals, 'delimiter', ' ', 'precision', 6,'-append');
    
end