clear all;
close all;

addpath('export_fig/');
numPts = 25;
%path = ['../results/Non-stationary/' num2str(numPts) '/'];

path = ['../results/WLake/' num2str(numPts) '/No-Corrections'];
%path = ['../results/WLake/200/MinConditioning/'];
listing = dir(path);

for i = 3:length(listing)
    filename = [path listing(i).name];
    
    filename = '/home/lewisli/cuda-workspace/MPSKriging/results/WLake/25/No-Corrections/MultiTI3_PenaltyOff_MinCond1_MaxCond50_Search25.dat'
    
    A = importdata(filename);
    estimate = A(:,3);
    variance = A(:,4);
    
    gridSize = [max(A(:,2))+1 max(A(:,1))+1];
    saveName = [listing(i).name num2str(numPts) 'Est.png']
    displayRealizations(estimate,gridSize, [listing(i).name '' ...
        num2str(numPts) ' Hard Data Pts Estimate']) ;
    % export_fig([path saveName],'-png', '-m2','-nocrop');
    figure;
    displayRealizations(variance,gridSize, [listing(i).name '' ...
        num2str(numPts) ' Hard Data Pts Variance']) ;
    
    break;
end
% filename = '../results/GaussianExample/25/TI5Results50';


%
% %
% figure
% displayRealizations(variance,gridSize, [filename ' Variance']) ;
% saveName = [filename 'Var.png' ];
% %export_fig(saveName,'-png', '-m2');
%close all;

%%break;


%%end



