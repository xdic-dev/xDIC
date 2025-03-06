
%%% compute_resolution_pointcloud 
% Take the closest point and compute its mean distance for all points

function varargout = compute_resolution_pointcloud(points,showgraph)
    numNeighbours = 1; 
    %create kdtree
    kdtreeobj = KDTreeSearcher(points,'distance','euclidean');
    %get nearest neighbours
    n = knnsearch(kdtreeobj,points,'k',(numNeighbours+1));
    %remove self
    n = n(:,2:end);
    %find absolute difference in position from neighbouring points
    p = abs(repmat(points(:,1:3),numNeighbours,1) - points(n(:),1:3));
    p = reshape(p, size(points,1),numNeighbours,3);
    %compute the mean distance from all neigbouring points 
    dist = squeeze(mean(p,2,'omitnan')); 
    %compute the norm of euclidian distance 
    dist_norm = sqrt(dist(:,1).^2+dist(:,2).^2+dist(:,3).^2); 
    
    %graph
    if showgraph
        newfig; 
        Nrow = 1; Ncol = 1;  
        ax = subplot_ax(Nrow,Ncol); hold(ax,'on'); 
        boxplot(dist_norm);
        title(ax,sprintf('Distribution of distance to neighborhood points\nmean = %1.3fmm',mean(dist_norm,'omitnan'))); 
        ax.XTickLabel = [];
        set(get(ax, 'XAxis'), 'Visible', 'off');
        ax.YLabel.String = 'distance [mm]'; 
        set(ax,'box','off'); 
        ax.XLim = [0.80,1.20]; 
        ax.YLim = [0,1];
    end
    
    if nargout >= 1
        varargout{1} = dist_norm; 
    end
end