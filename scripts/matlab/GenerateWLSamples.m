% GenerateWLSamples
% Sample hard data points from Walker Lake
%
% Author: Lewis Li
% email: lewisli@stanford.edu
% Website: Website: http://scrf.stanford.edu/
% Feb 2015; Last revision: 28-Feb-2015

% ------------- ESTIMATION PARAMETERS --------------
trueDEMPath = '../../data/DS-NonStationary/';
caseName = 'DS-NonStationary';
trueDEM = 'Reference';

% trueDEMPath = '../../data/WLake/';
% trueDEM = 'WalkerLake';

SYNCDATAPATH='/Volumes/Data/Data/UK-TI/hardData/';

[TruthDim TruthVals TruthNames] = loadObjects([trueDEMPath trueDEM]);
displayRealizations(TruthVals,TruthDim,TruthNames);

%%

numPts = [400 600];
numSamples = 50;

for i = 1:length(numPts)
    for j= 1:numSamples
        numHardPoints = numPts(i);
        HDatVals = sampleHardData(TruthVals(:,1),TruthDim, numHardPoints);
        HDatVals(:,1:3) = HDatVals(:,1:3) - 1;

        outputDir = [SYNCDATAPATH caseName '/' num2str(numPts(i))];
        outputFileName = [outputDir '/hData' num2str(j)]

        fid=fopen(outputFileName, 'w');
        txt = [ 'random samples\n4\nX\nY\nZ\nU-random-sampling\n'];
        S = sprintf(txt);
        fprintf(fid,'%s',S);
        
        dlmwrite([outputFileName], HDatVals, 'delimiter', ' ', 'precision', 6,'-append');
    end
end
