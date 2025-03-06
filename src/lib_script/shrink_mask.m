function [mask_out,B] = shrink_mask(mask, factor)
siz = size(mask); 
B = bwboundaries(imgradient(mask),'holes');
A = B{1}; 
xi = A(:,1);
yi = A(:,2);
j = boundary(xi,yi); 
xi = xi(j);
yi = yi(j);
p = polyshape(xi,yi);
[cx, cy] = centroid(p); %find the centroid
m = scale(p,factor,[cx cy]);
xi = m.Vertices(:,1);
yi = m.Vertices(:,2);
mask_out = logical(poly2mask(yi,xi,siz(1),siz(2)));
end