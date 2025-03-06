function data_cat = cat_imstruct_data(data_struct,field,varargin)
    % Parse input arguments
    p = inputParser;
    p.addRequired('data_cell');
    p.addRequired('field');
    p.addParameter('numfield',1);
    p.addParameter('idx',[]);
    
    % Parse input arguments
    p.parse(data_struct,field,varargin{:});
    numfield = p.Results.numfield;
    idx = p.Results.idx;
    
    if ~isempty(idx)
        Ntrial = length(find(idx == 1));
        mask = find(idx == 1); 
    else
        Ntrial = length(data_struct); 
        mask = 1:Ntrial;
    end
    
    data_cat = zeros(length(data_struct(mask(1)).(field)),Ntrial);
    for ii = 1:Ntrial
        data_cat(:,ii) = data_struct(mask(ii)).(field)(:,numfield);
    end
end