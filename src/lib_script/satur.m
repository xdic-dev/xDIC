function var = satur(var, varargin)
% Parse inputs
p = inputParser;
p.addParameter('level',200);
p.addParameter('method','high');

p.parse(varargin{:});
level = p.Results.level;
method = p.Results.method;

if length(level)~=1 
    Nvar = length(level);
else 
    Nvar = 1; 
end
for ii = 1:Nvar
    if strcmp(method,'high')
        var(var > level(ii)) = level(ii);
    elseif strcmp(method,'low')
        var(var < level(ii)) = level(ii);
    end
end
end