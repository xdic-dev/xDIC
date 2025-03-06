function data_cat = cat_imstruct_data(data_struct,field,varargin)
    % Parse input arguments
    p = inputParser;
    p.addRequired('data_struct');
    p.addRequired('field');
    p.addParameter('idx',[]);
    
    % Parse input arguments
    p.parse(data_struct,field,varargin{:});
    idx = p.Results.idx;
    
    if ~isempty(idx)
        Ntrial = length(find(idx == 1));
        mask = find(idx == 1); 
    else
        Ntrial = length(data_struct); 
        mask = 1:Ntrial;
    end
    siz = size(data_struct(mask(1)).(field));
    data_cat = zeros(siz(1),siz(2),Ntrial);
    for ii = 1:Ntrial
        data_cat(:,:,ii) = data_struct(mask(ii)).(field);
    end
end