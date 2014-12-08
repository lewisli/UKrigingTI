% loadObjects - Function to load SGEMS realizations into MATLAB workspace
% SGEMS stores realizations by default in binary format which is hard to
% work with, so instead we use .gslib
% Syntax:  [dim realVals realNames] = loadObjects(path)
%
% Inputs:
%    path - path of the sgems object
%
% Outputs:
%    dim - dimensions of the grid
%    obj - actual object in matrix format
%  names - names of realizations
%
% Example:
%    [dim val name] = loadObjects("data/U.sgems)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%
% Author: Lewis Li
% email: lewisli@stanford.edu
% Website: http://www.
% July 2013; Last revision: 03-July-2013

%------------- BEGIN CODE --------------
function [dim realVals realNames] = loadObjects(path)
addpath('../data/DS-NonStationary');
%addpath('data/Gaussian');
%addpath('data/Non-stationary');
fid = fopen(path,'r');

% We need to differentiate between full grid and sampled points, the former
% has a dimension associated with it

% Check first line to see if it contains
% First line is the grid name and dimensions
InputText=textscan(fid,'%s',1,'delimiter','\n');
gridNameStr = InputText{1,1};

expression = '((\w*)';
dimStr = regexp(gridNameStr{1,1},expression,'match');

% If we are loading full grid
if (length(dimStr) > 0)
    dimStr = strrep(dimStr, 'x', ' ');
    dimStr = strrep(dimStr, '(', '');
    dim = str2num(dimStr{1,1});
    
    % Next line is number of realizations on grid
    InputText=textscan(fid,'%s',1,'delimiter','\n');
    numReal = str2double(InputText{1,1});
    
    % Next lines are names of realizations
    InputText=textscan(fid,'%s',numReal,'delimiter','\n');
    realNames = InputText{1,1};
    
    % Read in matrix of data
    rawData = textscan(fid, '%f', 'delimiter', ' ');
    realVals = cell2mat(rawData(1,1));
    realVals = reshape(realVals,numReal,[])';
    
    %displayRealizations(realVals(:,1),dim,'W');
else
    % If we are loading samples then the first row indicates the name and
    % the second row indicates the number of properties we have sampled
    InputText=textscan(fid,'%s',1,'delimiter','\n');
    numReal = str2double(InputText{1,1});
    
    % Next lines are names of realizations
    InputText=textscan(fid,'%s',numReal,'delimiter','\n');
    realNames = InputText{1,1};
    
    % Read in matrix of data
    rawData = textscan(fid, '%f', 'delimiter', ' ');
    realVals = cell2mat(rawData(1,1));
    realVals = reshape(realVals,numReal,[])';
    
    % First 3 columns are 0 indexed
    realVals(:,1:3) = realVals(:,1:3) +1;
    
    % Get number of hard data points
    dim = length(realVals);
end


end