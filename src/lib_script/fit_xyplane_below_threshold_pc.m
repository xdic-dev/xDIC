function varargout = fit_xyplane_below_threshold_pc(points,threshold,varargin)
    
    %% Parse input arguments
    p = inputParser;
    p.addRequired('points');
    p.addRequired('threshold');
    p.addParameter('margin',0); 
    p.parse(points,threshold,varargin{:});
    margin = p.Results.margin;
    
    %mask
    below_logic = points(:,3) < threshold+margin;
    %indices
    below_idx = find(below_logic); 
    [temp,fun_planefit] = fit_xyplane_pc(points(below_logic,:),'margin',margin);
    %mask recomputed
    mask_in = false(size(points,1),1); 
    mask_in(below_idx(temp,:)) = true; 
    
    
    %% outputs arguments
    varargout{1} = mask_in; 
    varargout{2} = fun_planefit; 


end
