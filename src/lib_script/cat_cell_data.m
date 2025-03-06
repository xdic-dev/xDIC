function data_cat = cat_cell_data(data_cell,field,varargin)
    % Parse input arguments
    p = inputParser;
    p.addRequired('data_cell');
    p.addRequired('field');
    p.addParameter('numfield',1);
    p.addParameter('mask',[]);
    p.addParameter('idx',[]);
    
    % Parse input arguments
    p.parse(data_cell,field,varargin{:});
    numfield = p.Results.numfield;
    mask = p.Results.mask;
    idx = p.Results.idx;
    
    if ~isempty(mask)
        Ntrial = length(find(mask == 1));
        mask = find(mask == 1); 
    else
        Ntrial = length(data_cell);
        mask = 1:Ntrial;
    end
    if isempty(idx)
        Ndata = size((data_cell{mask(1)}.(field)),1); 
        idx = 1:Ndata;
    else
        Ndata = length(idx);
    end
    data_cat = zeros(Ndata,Ntrial);
    for ii = 1:Ntrial
        data_cat(idx,ii) = data_cell{mask(ii)}.(field)(idx,numfield);
    end
end