function varargout = map_pixel2subset(pos,spacing)
    varargout{1} = round((pos-1)/(spacing+1)+1);
end