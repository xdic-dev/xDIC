function varargout = map_pointcoordinate(pos, field)
%field
u = field{1};
v = field{2};
%position
x = pos(1);
y = pos(2);
varargout{1} = round([x+u(y,x),y+v(y,x)]);
end