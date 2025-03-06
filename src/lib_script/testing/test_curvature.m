
newfig('test measurement curvature'); 


clim_set = 0.50;
npoint = 150; 

subplot(1,2,1); hold on; 
iframe = 1; 
ptFC_cur = DIC3Dcombined.Points3D{iframe}; 
%Look at the distance between neighborhood points and detect
%outliers according to a multiple (threshold) of the variance.
[~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(ptFC_cur),...
    'NumNeighbors',15,...
    'Threshold',1);
temp = false(length(ptFC_cur),1);
temp(outlierIndices) = true;

x = ptFC_cur(:,1);
y = ptFC_cur(:,2);
z = ptFC_cur(:,3);
% figure; hold on; 
% scatter3(x(temp),y(temp),z(temp),'MarkerFaceColor','r','MarkerEdgeColor','none');
% scatter3(x(~temp),y(~temp),z(~temp),'MarkerFaceColor','b','MarkerEdgeColor','none');

logicNaN = isnan(x);
mask_outliers = temp|logicNaN;

xg = linspace(-15, max(x), npoint); % coordinates of grid
yg = linspace(-10, max(y), npoint);
[X, Y] = meshgrid(xg, yg); % X and Y grids
% Z = griddata(x(~mask_outliers), y(~mask_outliers), z(~mask_outliers), X, Y);
Z = griddata(x(~mask_outliers), y(~mask_outliers), z(~mask_outliers), X, Y);
% figure; surf(X, Y, Z);
% set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');
[K,H,P1,P2] = surfature(X,Y,Z);
surf(X,Y,Z,H,'facecolor','interp','EdgeColor','k');
set(gca,'clim',[-clim_set,clim_set])
set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000]);
colormap jet
colorbar

subplot(1,2,2); hold on; 
iframe = 150; 
ptFC_cur = DIC3Dcombined.Points3D{iframe}; 
%Look at the distance between neighborhood points and detect
%outliers according to a multiple (threshold) of the variance.
[~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(ptFC_cur),...
    'NumNeighbors',15,...
    'Threshold',1);
temp = false(length(ptFC_cur),1);
temp(outlierIndices) = true;

x = ptFC_cur(:,1);
y = ptFC_cur(:,2);
z = ptFC_cur(:,3);
% figure; hold on; 
% scatter3(x(temp),y(temp),z(temp),'MarkerFaceColor','r','MarkerEdgeColor','none');
% scatter3(x(~temp),y(~temp),z(~temp),'MarkerFaceColor','b','MarkerEdgeColor','none');

logicNaN = isnan(x);
mask_outliers = temp|logicNaN;

xg = linspace(-15, max(x), npoint); % coordinates of grid
yg = linspace(-10, max(y), npoint);
[X, Y] = meshgrid(xg, yg); % X and Y grids
% Z = griddata(x(~mask_outliers), y(~mask_outliers), z(~mask_outliers), X, Y);
Z = griddata(x(~mask_outliers), y(~mask_outliers), z(~mask_outliers), X, Y);
% figure; surf(X, Y, Z);
% set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');
[K,H,P1,P2] = surfature(X,Y,Z);
surf(X,Y,Z,H,'facecolor','interp','EdgeColor','k');
set(gca,'clim',[-clim_set,clim_set])
set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000]);
colormap jet
colorbar

%% 3D curvature - surface variance

clim_set = [0,0.025];%0.05; 

gridStep = 0.75;%0.75;
numNeighbours = 50; 

id = [1,150]; 
scatter_width = 15;%15;

h_fig = newfig(sprintf('Surface variance S9 %1.3f %1.3f',gridStep,numNeighbours)); 
set(h_fig,'Position',[0.0265    0.8202   19.5527   25.5058]);
ii = 0; 
for jj = 1:3
for iframe = id
ii = ii+1; 
subplot(3,length(id),ii); hold on; 
ptFC_cur = DIC3Dcombined.Points3D{iframe}; 
% ptFC_cur = DIC3DPPresults.Points3D{iframe}; 
%Look at the distance between neighborhood points and detect
%outliers according to a multiple (threshold) of the variance.
ptCloud = pointCloud(ptFC_cur); 
[~,inlierIndices1,outlierIndices1] = pcdenoise(ptCloud,...
    'NumNeighbors',15,...
    'Threshold',3);
[~,inlierIndices2,outlierIndices2] = pcdenoise(pointCloud(ptCloud.Location(inlierIndices1,:)),...
    'NumNeighbors',15,...
    'Threshold',3);
temp = true(length(ptFC_cur),1);
temp(inlierIndices1(inlierIndices2)) = false;

logicNaN = isnan(ptFC_cur(:,1));
mask_outliers = logicNaN|temp; 
points = ptFC_cur(~mask_outliers,:);
% Points within each grid box are merged by averaging their locations
ptCloud = pointCloud(points); 
ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);
pointsds = ptCloudA.Location;
% stepSize = floor(ptCloud.Count/ptCloudA.Count);
% indices = 1:stepSize:ptCloud.Count;
% pointsds = points(indices,:);

% 
temp = pointCloud(pointsds); 
temp = pcdownsample(temp,'gridAverage',gridStep);
pointsds = temp.Location;

temp = pointCloud(pointsds); 
temp = pcdownsample(temp,'gridAverage',gridStep);
pointsds = temp.Location;

temp = pointCloud(pointsds); 
temp = pcdownsample(temp,'gridAverage',gridStep);
pointsds = temp.Location;

%find the normals and curvature
[normals,curvature] = findPointNormals(pointsds,numNeighbours,[0,0,100],true);
scatter3(points(:,1),points(:,2),points(:,3),1,'filled'); 
scatter3(pointsds(:,1),pointsds(:,2),pointsds(:,3),scatter_width,curvature,'filled'); 

title_set = sprintf('#%d',iframe); 
view_plot = [0 -90];

if jj == 1
    title_set = sprintf('Y-Z #%d',iframe); 
    view_plot = [90 0];
elseif jj == 2
    title_set = sprintf('X-Z #%d',iframe); 
    view_plot = [0 0];
elseif jj == 3
    title_set = sprintf('X-Y #%d',iframe); 
    view_plot = [0 -90];
end
set(gca,'clim',clim_set)
xlim(gca,[-15,15]); 
ylim(gca,[-10,25]); 
zlim(gca,[-2,15]); 
title(gca,title_set); 
set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', view_plot);
colormap summer
colorbar
end
end
%% 2D curvature
clim_set = [0,0.01];%0.05; 

gridStep = 1;%0.75;
numNeighbours = 8; 

id = [30]; 
scatter_width = 50;%15;

h_fig = newfig(sprintf('Measurement of curvature 12 %1.3f %1.3f',gridStep,numNeighbours)); 
set(h_fig,'Position',[0.0265    0.8202   19.5527   25.5058]);
ii = 0; 
for axis_target = [1]%,2,3]
    for iframe = id
        ii = ii+1; 
        subplot(1,length(id),ii); hold on; 
        % ptFC_cur = DIC3Dcombined.Points3D{iframe}; 
        ptFC_cur = DIC3DPPresults.Points3D{iframe}; 
        %Look at the distance between neighborhood points and detect
        %outliers according to a multiple (threshold) of the variance.
        ptCloud = pointCloud(ptFC_cur); 
        [~,inlierIndices1,outlierIndices1] = pcdenoise(ptCloud,...
            'NumNeighbors',15,...
            'Threshold',3);
        [~,inlierIndices2,outlierIndices2] = pcdenoise(pointCloud(ptCloud.Location(inlierIndices1,:)),...
            'NumNeighbors',15,...
            'Threshold',3);
        temp = true(length(ptFC_cur),1);
        temp(inlierIndices1(inlierIndices2)) = false;

        logicNaN = isnan(ptFC_cur(:,1));
        mask_outliers = logicNaN|temp; 
        points = ptFC_cur(~mask_outliers,:);
        % Points within each grid box are merged by averaging their locations
        ptCloud = pointCloud(points); 
        ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);
        pointsds = ptCloudA.Location;
        % stepSize = floor(ptCloud.Count/ptCloudA.Count);
        % indices = 1:stepSize:ptCloud.Count;
        % pointsds = points(indices,:);

        other_axis = 1:3; other_axis(axis_target) = []; 
        % range_point = max(pointsds(:,axis_target))-min(pointsds(:,axis_target));
        % idx_mask = find(pointsds(:,axis_target)<mean(pointsds(:,axis_target))-0.05*range_point&pointsds(:,axis_target)>mean(pointsds(:,axis_target))-0.35*range_point); 
        pointsds(:,axis_target) = mean(pointsds(:,axis_target)); 
        % pointsds = pointsds(idx_mask,:);

        temp = pointCloud(pointsds); 
        temp = pcdownsample(temp,'gridAverage',gridStep);
        pointsds = temp.Location;

        k = boundary(pointsds(:,other_axis(1)),pointsds(:,other_axis(2)),0.1); 
        pointsds = pointsds(k,:);
        [pointsds,dudt,fofthandle] = interparc(50,pointsds(:,1),pointsds(:,2),pointsds(:,3),'linear');


        %find the normals and curvature
        [normals,curvature] = findPointNormals(pointsds,numNeighbours,[0,0,100],true);
        angle = atan(abs(normals(:,other_axis(1)))./abs(normals(:,other_axis(2))))*180/pi; 
        scatter3(points(:,1),points(:,2),points(:,3),1,'filled'); 
        scatter3(pointsds(:,1),pointsds(:,2),pointsds(:,3),scatter_width,curvature,'filled'); 

        if axis_target == 1
            title_set = sprintf('Y-Z #%d',iframe); 
            view_plot = [90 0];
        elseif axis_target == 2
            title_set = sprintf('X-Z #%d',iframe); 
            view_plot = [0 0];
        elseif axis_target == 3
            title_set = sprintf('X-Y #%d',iframe); 
            view_plot = [0 -90];
        end
        set(gca,'clim',clim_set)
        xlim(gca,[-15,15]); 
        ylim(gca,[-10,25]); 
        zlim(gca,[-2,15]); 
        title(gca,title_set); 
        set(gca,'DataAspectRatio',[1 1 1],'Color','None','View',view_plot);
        colormap(parula);
        colorbar
    end
end

curv_lim_set = 0.0005;
zlim_set = 0; 
xvec = -10:25; 
mask_curv = curvature < curv_lim_set; 
z = pointsds(:,3); 
size(z)
size(mask_curv)
mask1 = mask_curv & z>zlim_set; 
mask2 = mask_curv & z<zlim_set; 

%mask 1
x = pointsds(mask1,2); y = pointsds(mask1,3);
c1 = polyfit(x,y,1); y_est1 = polyval(c1,xvec);
%mask 2
x = pointsds(mask2,2); y = pointsds(mask2,3);
c2 = polyfit(x,y,1); y_est2 = polyval(c2,xvec);
angle = (atan(c1(1))-atan(c2(1)))*180/pi
plot3(gca,zeros(length(xvec),1),xvec,y_est1,'k');
plot3(gca,zeros(length(xvec),1),xvec,y_est2,'k');



%%
iframe = 1; 
gridStep = 0.5; 

ptFC_cur = DIC3Dcombined.Points3D{iframe}; 
%Look at the distance between neighborhood points and detect
%outliers according to a multiple (threshold) of the variance.
ptCloud = pointCloud(ptFC_cur); 
[~,inlierIndices1,outlierIndices1] = pcdenoise(ptCloud,...
    'NumNeighbors',15,...
    'Threshold',3);
[~,inlierIndices2,outlierIndices2] = pcdenoise(pointCloud(ptCloud.Location(inlierIndices1,:)),...
    'NumNeighbors',15,...
    'Threshold',3);
temp = true(length(ptFC_cur),1);
temp(inlierIndices1(inlierIndices2)) = false;

logicNaN = isnan(ptFC_cur(:,1));
mask_outliers = logicNaN|temp; 
points = ptFC_cur(~mask_outliers,:);

ptCloudA = pcdownsample(pointCloud(points),'gridAverage',gridStep);
stepSize = floor(ptCloud.Count/ptCloudA.Count);
indices = 1:stepSize:ptCloud.Count;
ptCloudB = select(ptCloud,indices);
figure
subplot(1,2,1); pcshow(ptCloudA);
subplot(1,2,2); pcshow(ptCloudB);

% step_length = 0.75; 
% xq = min(x):step_length:max(x); % coordinates of grid
% yq = min(y):step_length:max(y); % coordinates of grid
% zq = min(z):step_length:max(z); % coordinates of grid
% vq = griddata(x,y,z,v,xq,yq,zq);

%plot normals and colour the surface by the curvature
% hold off;
% surf(x,y,z,reshape(curvature,49,49));
% hold on;
% quiver3(points(:,1),points(:,2),points(:,3),...
%     normals(:,1),normals(:,2),normals(:,3),'r');
% axis equal;
%%
ptCloud = pointCloud([x,y,z]);
gridStep = 1 ;
ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);

figure; 
subplot(2,1,1); pcshow(ptCloud)
subplot(2,1,2); pcshow(ptCloudA)
%%
[X,Y] = meshgrid(min(x):max(x));
%%
% [X,Y,Z] = peaks; 

figure; 

%Qucik Program to demo the use of findPointNormals
%generate a set of 3d points
x = repmat(1:49,49,1);
y = x';
z = peaks;
points = [x(:),y(:),z(:)];
%find the normals and curvature
[normals,curvature] = findPointNormals(points,[],[0,0,10],true);
%plot normals and colour the surface by the curvature
hold off;
surf(x,y,z,reshape(curvature,49,49));
hold on;
quiver3(points(:,1),points(:,2),points(:,3),...
    normals(:,1),normals(:,2),normals(:,3),'r');
axis equal;