function varargout = smooth_space(FC,P,optStruct,varargin)
%% Parse input arguments
p = inputParser;
p.addRequired('FC');
p.addRequired('P');
p.addRequired('optStruct');
p.addParameter('idx_nearest_neighbours_set',[]);
p.addParameter('weight_nearest_neighbours_set',[]);
p.parse(FC,P,optStruct,varargin{:});

idx_nearest_neighbours_set = p.Results.idx_nearest_neighbours_set;
weight_nearest_neighbours_set = p.Results.weight_nearest_neighbours_set;
numNeighbours = optStruct.smoothPar.n; % SPACE FILTERING PARAM 1
sigma = optStruct.smoothPar.sigma; % SPACE FILTERING PARAM 2

%% setup
%ensure inputs of correct type
P = double(P);
if isempty(idx_nearest_neighbours_set)
    %create kdtree : knnsearch and rangesearch can be used 
    kdtreeobj = KDTreeSearcher(P,'distance','euclidean');
    %get k-nearest neighbors : non parametric supervised learning method 
    n = knnsearch(kdtreeobj,P,'k',(numNeighbours+1));
    
    %find difference in position from neighbouring points
    p = repmat(P(:,1:3),(numNeighbours+1),1) - P(n(:),1:3);
    p = mean(reshape(p, size(P,1),(numNeighbours+1),3),3,'omitnan');
    
    % weight from distance
    dist = compute_resolution_pointcloud(P,0); dist = mean(dist,'omitnan');
    W = exp(-((p.*(1/dist)).^2./sigma^2));
    W = W./sum(W,2);
else
    n = idx_nearest_neighbours_set;
    W = weight_nearest_neighbours_set;
end


%% Find face value of neighbours 
FCsmooth = FC(n(:));
FCsmooth = reshape(FCsmooth, size(P,1),(numNeighbours+1));

FCsmooth = sum(FCsmooth .* W,2,'omitnan'); 
%% output
if nargout >= 1
    varargout{1} = FCsmooth;
end
if nargout >= 2
    varargout{2} = n;
end
if nargout >= 3
    varargout{3} = W;
end
end