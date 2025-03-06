function varargout = find_curv_threshold(points,varargin)

% idx_target_curv is an indice that allows to make 2D curvature measurement
% using the same method. If you have 3d point cloud and squeesh a dimension to 
% get the curvature in one axes. 

%% Parse input arguments
p = inputParser;
p.addRequired('points'); 
p.addParameter('numNeighbours',9);
p.addParameter('threshold_percent',[]);
p.addParameter('threshold_absolute',[]);
p.addParameter('idx_target_curv',1);
p.parse(points,varargin{:});
numNeighbours = p.Results.numNeighbours;
threshold_percent = p.Results.threshold_percent;
threshold_absolute = p.Results.threshold_absolute;
idx_target_curv = p.Results.idx_target_curv;

%%
nPoints = size(points,1); 
logicNaN = isnan(points(:,1));
idx_valid = find(~logicNaN);
[~,curvature] = findPointNormals(points(idx_valid,:),numNeighbours,[0,0,100],true);
curvature = curvature(:,idx_target_curv);
curvature_all = NaN(nPoints,1); 
curvature_all(idx_valid) = curvature; 
if ~isempty(threshold_percent)
    theshold = prctile(curvature,threshold_percent);
else
    theshold = threshold_absolute;
end
idx_curv = find(curvature > theshold);
mask_curv = false(nPoints,1);
mask_curv(idx_valid(idx_curv)) = true;

%% Output arguments 
if nargout >= 1
    varargout{1} = mask_curv;
end
if nargout >= 2
    varargout{2} = theshold;
end
if nargout >= 3
    varargout{3} = curvature_all;
end
end