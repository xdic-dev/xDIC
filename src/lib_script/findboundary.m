
function varargout = findboundary(p_raw,varargin)
    h_parser = inputParser;
    %this parameter allows to bound the extreme points more or less compact
    %A very compact boundary will be created with a value 1, meaning that 
    %the boundary will very likely be uneven/irregular.
    h_parser.addParameter('connex_factor',0.2);
    h_parser.parse(varargin{:});
    connex_factor = h_parser.Results.connex_factor;
    
    %filtering NaN
    idx = ~isnan(p_raw(:,1));
    p = p_raw(idx,:); %remove NaN points 
%     p(p(:,1) == 0,:) = []; %remove point at the origin
    
    %find boundary 
    x=p(:,1); y=p(:,2); 
    k = boundary(x,y,connex_factor);
    if size(p,2)==3 
        %if there is a z-coordinates
        z=p(:,3); 
        ca = [x(k),y(k),z(k)];
    else
        ca = [x(k),y(k)];
    end
    
    %retreive border indices
    mask_border = false(size(p_raw,1),1); 
    mask_border(idx(k)) = true; 
    
    %outputs
    if nargout >= 1 
        varargout{1} = ca; 
    end
    if nargout >= 2
       varargout{2} = mask_border;
    end
    
end