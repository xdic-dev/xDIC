p = DIC3DPPresults.Points3D{1}; 

x = p(:,1); 
y = p(:,2); 
z = p(:,3); 
[z0,idx] = min(z,[],1,'omitnan');
x0 = x(idx);
y0 = y(idx);

nanLogic = isnan(x); 
[~,~,V]=svd([zeros(length(x(~nanLogic)),1)-x0, y(~nanLogic)-y0, z(~nanLogic)-z0],0);
 a=V(1,end);
 b=V(2,end);
 c=V(3,end);
 d=a*x0+b*y0+c*z0;

 [Y,Z] = meshgrid(-10:20,0:15); 
figure; 
pcshow(pointCloud(p)); hold on; 
surf((-b/a)*Y+(-c/a)*Z+(d/a),Y,Z); hold on; 
plot3(x0,y0,z0,'r.','MarkerSize',20); 
