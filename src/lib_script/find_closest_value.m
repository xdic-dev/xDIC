function idx = find_closest_value(f,val)

val_format = zeros(1,length(val));
f_format = zeros(length(f),1);
val_format(:) = val; f_format(:) = f; 

p = abs(f_format - val_format);
[idx,~] = find(min(p) == p);

end