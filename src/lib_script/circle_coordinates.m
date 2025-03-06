function [x,y] = circle_coordinates(center,radius)
    angle = linspace(0,2*pi,1e3); 
    xi = center(1)+radius.*cos(angle); 
    yi = center(2)+radius.*sin(angle);
    xy = unique([ceil(xi'),ceil(yi')],'rows'); 
    x = xy(:,1); 
    y = xy(:,2); 
end