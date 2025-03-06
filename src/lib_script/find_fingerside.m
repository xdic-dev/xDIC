


function varargout = find_fingerside(PCfirst,PClast,CAlast,Nregion)

%% fit ellipse from last contact contour
[CAfit,centerEl] = myfitellipse(CAlast(:,1),CAlast(:,2)); 
%
CAfit_1 = (CAfit - centerEl)*(1/3)+centerEl; 
CAfit_2 = (CAfit - centerEl)*(2/3)+centerEl; 
CAfit_3 = CAfit; 

% indices of points inside the contour segmentation 
idx_IN1 = inpolygon(PClast(:,1),PClast(:,2),CAfit_1(:,1),CAfit_1(:,2));
idx_IN2 = inpolygon(PClast(:,1),PClast(:,2),CAfit_2(:,1),CAfit_2(:,2))&~idx_IN1;
idx_IN3 = inpolygon(PClast(:,1),PClast(:,2),CAfit_3(:,1),CAfit_3(:,2))&~idx_IN2&~idx_IN1;

zplane_mean = mean(CAlast(:,3)); 
% zmargin = 3; %mm of margin above the plane of contact
% idx_transversePlane = PClast(:,3)<zplane_mean+zmargin;%CAfit(max(PClast(:,3)))*0.33; 
idx_transversePlane2 = PClast(:,3)<max(PClast(:,3))*0.75; 
CAfit = [CAfit,repelem(zplane_mean,size(CAfit,1),1)];

%% Find center of fingertip contact 
% mean_displ = mean(PClast,1,'omitnan') - mean(PCfirst,1,'omitnan'); 
center = [mean(CAfit(:,1)),mean(CAfit(:,2)),mean(CAfit(:,3))]; 
% center = center + mean_displ; 

%% find symmetrical plane with the first 3d point cloud
% mask = true(size(PCfirst,1)); %~idx_IN; 
% DM = [PCfirst(:,2),PCfirst(:,3),ones(length(PCfirst(:,1)),1)]; % Design Matrix
% param = DM\PCfirst(:,1); % Estimate Parameters
% B = param; 
% fit_sagittalPlane = @(y,z) B(1)*y + B(3); %+ B(2)*z

% Translation to put the center at origin
trans = -[center(1) center(2) center(3)];
tform_translation = myrigidtform3d([0 0 0],trans);

%Point cloud replaced
pc_centered = pctransform(pointCloud(PClast),tform_translation);
%Contour replaced
calast_centered = pctransform(pointCloud(CAlast),tform_translation);

%loc
pc_centered_loc = pc_centered.Location;
calast_centered_loc = calast_centered.Location;

nanLogic = isnan(pc_centered_loc(:,1)); 
x0 = 0; y0 = 0; z0 = 0; 
[~,~,V]=svd([zeros(length(pc_centered_loc(~nanLogic,1)),1)-x0, pc_centered_loc(~nanLogic,2)-y0, pc_centered_loc(~nanLogic,3)-z0],0);
a=V(1,end); b=V(2,end); c=V(3,end);
d=a*x0+b*y0+c*z0;
%if we get rid of the z component in the plane 
fit_sagittalPlane = @(y,z) (-b/a)*y+(d/a);
%REMARK : We remove the z component of the plane fit, the height of the reconstructed point 
% can differ from a pair to the other 
% --> (-c/a)*z removed

% figure;
% ax=gca; hold(ax,'on');
% scatter3(ax,pc_centered_loc(:,1),pc_centered_loc(:,2),pc_centered_loc(:,3),10,'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
% scatter3(ax,fit_sagittalPlane(pc_centered_loc(:,2),pc_centered_loc(:,3)),pc_centered_loc(:,2),pc_centered_loc(:,3),20,'MarkerEdgeColor','b','MarkerFaceColor','b');
% scatter3(ax,0,0,0,20,'MarkerEdgeColor','r','MarkerFaceColor','r');
% scatter3(ax,calast_centered_loc(:,1),calast_centered_loc(:,2),calast_centered_loc(:,3),20,'MarkerEdgeColor','k','MarkerFaceColor','k');
% axis equal

%% Region segmentation
% Front AREA 
angle = [+90 0 0];
tform = myrigidtform3d(angle,[0 0 0]);
pCtransform = pctransform(pc_centered,tform);
idx_frontalPlane = pCtransform.Location(:,1) < fit_sagittalPlane(pCtransform.Location(:,2),pCtransform.Location(:,3)); 

alpha = 180/Nregion; %front area arc angle 

% In practice, we rotate the fingertip... so we rotate by the opposite
% amount.
idx_reg = false(length(pc_centered_loc),Nregion);
for ireg = 1:Nregion
    angular_rotation_finger = -90+alpha*ireg; 
    tform = myrigidtform3d([angular_rotation_finger 0 0],[0 0 0]);
    pCtransform = pctransform(pc_centered,tform);
    idx_in_left_plane = pCtransform.Location(:,1) < fit_sagittalPlane(pCtransform.Location(:,2),pCtransform.Location(:,3)); 
    if ireg == 1
        idx_reg(:,ireg) = idx_in_left_plane&(idx_frontalPlane);
    else
        idx_reg(:,ireg) = idx_in_left_plane&(~idx_plane_before);
    end
    idx_plane_before = idx_in_left_plane; 
end

% %%
% angle = [+45 0 0];
% tform = myrigidtform3d(angle,[0 0 0]);
% pCtransform = pctransform(PCfirst_centered,tform);
% idx_plane1 = pCtransform.Location(:,1) <= sym_planefit(pCtransform.Location(:,2),pCtransform.Location(:,3)); 
% 
% angle = [45+90 0 0];
% tform = myrigidtform3d(angle,[0 0 0]);
% pCtransform = pctransform(PCfirst_centered,tform);
% idx_plane2 = pCtransform.Location(:,1) <= sym_planefit(pCtransform.Location(:,2),pCtransform.Location(:,3)); 
% %%
% angle = [+45 0 0];
% tform = myrigidtform3d(angle,[0 0 0]);
% pCtransform = pctransform(PCfirst_centered,tform);
% idx_plane1 = pCtransform.Location(:,1) > sym_planefit(pCtransform.Location(:,2),pCtransform.Location(:,3)); 
% 
% angle = [45+90 0 0];
% tform = myrigidtform3d(angle,[0 0 0]);
% pCtransform = pctransform(PCfirst_centered,tform);
% idx_plane2 = pCtransform.Location(:,1) < sym_planefit(pCtransform.Location(:,2),pCtransform.Location(:,3)); 
% %%
% C = colororder();
% colorR1 = C(1,:); %rgb('LightCoral');
% colorR2 = C(2,:); %rgb('CadetBlue');
% colorR3 = C(3,:); %rgb('PaleGreen');
% 
% figure;
% subplot(3,1,1); ax=gca; hold(ax,'on');
% scatter3(ax,PCfirst(idx_R1,1),PCfirst(idx_R1,2),PCfirst(idx_R1,3),10,'MarkerEdgeColor',colorR1,'MarkerFaceColor',colorR1);
% subplot(3,1,2); ax=gca; hold(ax,'on');
% scatter3(ax,PCfirst(idx_R2,1),PCfirst(idx_R2,2),PCfirst(idx_R2,3),10,'MarkerEdgeColor',colorR2,'MarkerFaceColor',colorR2);
% subplot(3,1,3); ax=gca; hold(ax,'on');
% scatter3(ax,PCfirst(idx_R3,1),PCfirst(idx_R3,2),PCfirst(idx_R3,3),10,'MarkerEdgeColor',colorR3,'MarkerFaceColor',colorR3);
% 
% plot3(ax,CAlast(:,1),CAlast(:,2),CAlast(:,3),'k');
% scatter3(ax,center(1),center(2),center(3)-2,100,'k','filled');
% axis equal
% set(ax,'Color','None');
% set(ax,'visible','off');
% set(ax,'View', [0 -90.0000]);

% set(gcf,'Color','None'); 
%% output 
if nargout >= 1
    varargout{1} = idx_reg; 
end
if nargout >= 2
    varargout{2} = idx_IN1; 
end
if nargout >= 3
    varargout{3} = idx_IN2; 
end
if nargout >= 4
    varargout{4} = idx_IN3; 
end
if nargout >= 5
    varargout{5} = ~idx_transversePlane2&~idx_frontalPlane; 
end
end

