function data_cat = cat_smp_data(data_robot,field,varargin)
    % Parse inputs
    [numfield,~] = parseargpair(varargin,'numfield',1);
    [idx,~] = parseargpair(varargin,'idx',[]);
    if ~isempty(idx)
        Ntrial = length(find(idx == 1));
    else
        Ntrial = length(data_robot); 
    end
    data_cat = zeros(length(data_robot{mask(1)}.time),Ntrial);
    for i = 1:Ntrial
        data_cat(:,i) = data_robot{i}.(field)(:,numfield);
    end
end