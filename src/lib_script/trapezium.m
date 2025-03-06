function coor = trapezium(a,b,h) 
% h = 1 ;  % height 
% a = 2 ;  % top side
% b = 4 ;   % base 
%%Frame vertices
A = [0 0] ;
B = [+a 0] ;
C = [+b h] ;
D = [0 h] ;  
coor = [A ; B; C; D] ;  
% patch(coor(:,1), coor(:,2),'r')
end