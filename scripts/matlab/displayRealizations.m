% displayRealizations - Function to plot realizations
%
% Syntax:  [h] = displayRealizations(dat, dim, name)
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

function [h] = displayRealizations(dat, dim, name)
maxVal = max(dat);
minVal = min(dat);
realization = reshape(dat,dim);
imagesc((rot90(realization)))
axis image;
h = gca;
colorbar;
set(h, 'XTick', []);
set(h, 'YTick', []);
t = title(name, 'FontSize', 34,'FontName', 'Helvetica-Narrow');
end

