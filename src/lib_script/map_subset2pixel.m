
function varargout = map_subset2pixel(pos,spacing)
    varargout{1} = (pos-1)*(spacing+1)+1;
end