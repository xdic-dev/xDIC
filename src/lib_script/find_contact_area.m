
function varargout = find_contact_area(PCnow,PClast,mask_curv_now,varargin)
    %% Parse input arguments
    p = inputParser;
    p.addRequired('PCnow');
    p.addRequired('PClast');
    p.addRequired('mask_curv_now');
    p.addParameter('zthreshold_lastmap',[]);
    p.addParameter('DispZnow',[]);
    p.addParameter('min_point_inside',3);
    p.addParameter('margin_around_plane',2);
    p.addParameter('numNeighboursCurv',20);
    p.addParameter('threshold_absolute_set',0.001);


    p.parse(PCnow,PClast,mask_curv_now,varargin{:});
    zthreshold_lastmap = p.Results.zthreshold_lastmap;
    DispZnow = p.Results.DispZnow;
    
    %allowable margin around fitted plane
    margin_around_plane = p.Results.margin_around_plane;      
    %Minimum number of point below threshold to be considered inside the
    %contact area
    min_point_inside = p.Results.min_point_inside; 
    
    %maximum possible indentation of the plateform (not completely rigid). 
    %if fitted plane is over this measure, all points are considered to be
    %outside the contact (or not already inside). 
    margin_platerform_inflection = 0.2; %[mm]
   
    %% Test to use DispZrate as proxi in the very start of loading 
    
    %% fit plane with points below threshold
    high_bound_threshold = zthreshold_lastmap+margin_platerform_inflection; 
        
    masknow_below = PCnow(:,3) < high_bound_threshold;
    mask_in = false(size(PCnow,1),1);  
    PCnow_below_id = false(size(PCnow,1),1);  
        
    if sum(masknow_below) > min_point_inside
        %first guess of points inside using plane fitting 
        temp = PCnow;
%         temp(PClast(:,2)<mean(PClast(:,2)),:) = NaN; %suppress a part of the fingertip 
        [PCnow_below_id,fun_planefit] = fit_xyplane_below_threshold_pc(temp,high_bound_threshold,'margin',margin_around_plane);

        % add curvature proxy 
        mask_in = PCnow_below_id&~mask_curv_now; %Attention : PCnow_below_id replaced
        
        %find contour
        if sum(mask_in) > min_point_inside
            %Look at the distance between neighborhood points and detect
            %outliers according to a multiple (threshold) of the variance.
            [~,inlierIndices,~] = pcdenoise(pointCloud(PClast(mask_in,:)),...
                'NumNeighbors',20,...
                'Threshold',3);
            temp = false(length(inlierIndices),1);
            temp(inlierIndices) = true;
            idx_in = find(mask_in);
            mask_in = false(size(mask_in)); 
            mask_in(idx_in(temp)) = true; 
        end
       
    end

    %% Parse output arguments
    if nargout >= 1
        varargout{1} = mask_in; 
    end
    if nargout >= 2
        varargout{2} = PCnow_below_id; 
    end
end

%     [X,Y] = meshgrid(min(PClast(:,1)):2:max(PClast(:,1)),min(PClast(:,2)):2:max(PClast(:,2)));
% 
%     temp1 = PClast_mask_in; 
%     Z1 = fun_planefit(X,Y); 
%     
%     [PClast_mask_in2,fun_planefit2] = fit_xyplane_below_threshold_pc(PClast,guess_threshold,'margin',margin_around_plane);
%     Z2 = fun_planefit2(X,Y); 
%     
%     [PClast_mask_in3,fun_planefit3] = fit_xyplane_below_threshold_pc(PClast,guess_threshold,'margin',2*margin_around_plane);
%     Z3 = fun_planefit3(X,Y); 
    
%     figure; 
%     subplot(2,3,1); hold on;     
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     scatter3(PClast(temp1,1),PClast(temp1,2),PClast(temp1,3),1,'g');
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     title('+0.0mm'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [90 0]);
%     
%     subplot(2,3,2); hold on;     
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     scatter3(PClast(PClast_mask_in2,1),PClast(PClast_mask_in2,2),PClast(PClast_mask_in2,3),1,'b');
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     surf(X,Y,Z2,'FaceAlpha',0.1); 
%     title('+0.1mm'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [90 0]);
%     
%     subplot(2,3,3); hold on;     
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     scatter3(PClast(PClast_mask_in3,1),PClast(PClast_mask_in3,2),PClast(PClast_mask_in3,3),1,'r');
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     surf(X,Y,Z3,'FaceAlpha',0.1);
%     title('+0.2mm'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [90 0]);
%     
%     subplot(2,3,3+1); hold on;     
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     scatter3(PClast(temp1,1),PClast(temp1,2),PClast(temp1,3),1,'g');
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     title('+0.0mm'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     subplot(2,3,3+2); hold on;     
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     scatter3(PClast(PClast_mask_in2,1),PClast(PClast_mask_in2,2),PClast(PClast_mask_in2,3),1,'b');
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     surf(X,Y,Z2,'FaceAlpha',0.1); 
%     title('+0.1mm'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     subplot(2,3,3+3); hold on;     
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     scatter3(PClast(PClast_mask_in3,1),PClast(PClast_mask_in3,2),PClast(PClast_mask_in3,3),1,'r');
%     surf(X,Y,Z1,'FaceAlpha',0.1);
%     surf(X,Y,Z3,'FaceAlpha',0.1);
%     title('+0.2mm'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
    
%     dist = compute_resolution_pointcloud(PClast,0); dist = mean(dist,'omitnan');
% %     compute_resolution_pointcloud(PCnow,1);
%     
%     figure; 
%     subplot(2,4,1); hold on; title('Point inside'); 
% %     surf(X,Y,Z1,'FaceAlpha',0.2); 
% %     surf(X,Y,Z2,'FaceAlpha',0.2); 
% %     surf(X,Y,Z3,'FaceAlpha',0.2);
%     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,[1 1 1].*0.1);
%     hp{3} = scatter3(PClast(PClast_mask_in3,1),PClast(PClast_mask_in3,2),PClast(PClast_mask_in3,3),3,'r','filled');
%     hp{2} = scatter3(PClast(PClast_mask_in2,1),PClast(PClast_mask_in2,2),PClast(PClast_mask_in2,3),3,'b','filled');
%     hp{1} = scatter3(PClast(temp1,1),PClast(temp1,2),PClast(temp1,3),3,'g','filled');
% %     plot3(PCnow(PCnow_in_border_id,1),PCnow(PCnow_in_border_id,2),PCnow(PCnow_in_border_id,3),'r');
% %     plot3(CAnow(:,1),CAnow(:,2),CAnow(:,3),'b','LineWidth',2);
%     h_leg = legend([hp{:}],{'+0mm','+1mm','+2mm'});  
%     h_leg.Box = 'off';
% %     h_leg.Title.String = 'Point inside'; 
% %     h_leg.Orientation = 'horizontal';
% %     set(h_leg,'FontSize',12); 
% 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     
%     ptile_curv = 50; 
%     %CURVATURE FROM ORIGINAL 
%     numNeighbours = 9;
%     scatter_width = 2;%15;
%     gridStep = 0.75; 
%     logicNaN = isnan(PClast(:,1));
%     idx1 = find(~logicNaN);
%     mask_below_threshold = PClast(idx1,3) < -1; 
%     [normals,curvature] = findPointNormals(PClast(idx1,:),numNeighbours,[0,0,100],true);
%     idx2 = find(curvature > prctile(curvature,ptile_curv)); 
%     mask_curv = false(size(PClast,1),1); 
%     mask_curv(idx1(idx2)) = true; 
%     clim_set = [0,prctile(curvature,ptile_curv)];
%     idx3 = find(PClast_mask_in2);    
%     
%     %Look at the distance between neighborhood points and detect
%     %outliers according to a multiple (threshold) of the variance.
%     [~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(PClast(PClast_mask_in3&~mask_curv,:)),...
%         'NumNeighbors',20,...
%         'Threshold',3);
%     temp = false(length(inlierIndices),1);
%     temp(inlierIndices) = true;
%     idx_in = find(PClast_mask_in3&~mask_curv);
%     PClast_mask_in3 = false(size(PClast_mask_in3)); 
%     PClast_mask_in3(idx_in(temp)) = true; 
%     
%     ca = findboundary(PClast(PClast_mask_in3&~mask_curv,:));
%     
%     % fit ellipse from last contact contour
%     [~,centerEl] = myfitellipse(ca(:,1),ca(:,2)); 
%     maskfront = ca(:,2)<centerEl(2); 
%     [cafit,~] = myfitellipse(ca(maskfront,1),ca(maskfront,2)); 
%     cafit = [cafit,repelem(mean(ca(:,3)),size(cafit,1),1)]; 
%     
%     
%     subplot(2,4,2); hold on; title(sprintf('Curvature\n from original pt'));
%     scatter3(PClast(idx1,1),PClast(idx1,2),PClast(idx1,3),scatter_width,curvature,'filled'); 
% %     scatter3(PClast(idx1(idx2),1),PClast(idx1(idx2),2),PClast(idx1(idx2),3),scatter_width,'r','filled'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     subplot(2,4,2+4); hold on; title(sprintf('Segmentation'));
%     scatter3(PClast(idx1,1),PClast(idx1,2),PClast(idx1,3),scatter_width,curvature,'filled'); 
%     hp1 = scatter3(PClast(PClast_mask_in3&~mask_curv,1),PClast(PClast_mask_in3&~mask_curv,2),PClast(PClast_mask_in3&~mask_curv,3),scatter_width,'b','filled'); 
%     hp2 = scatter3(PClast(mask_curv,1),PClast(mask_curv,2),PClast(mask_curv,3),scatter_width,'r','filled'); 
%     plot3(ca(:,1),ca(:,2),ca(:,3),'k-','LineWidth',3); 
% %     plot3(cafit(:,1),cafit(:,2),cafit(:,3),'k-','LineWidth',2,'Color','y'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     legend([hp1,hp2],{'in',sprintf('curv>%d%s',ptile_curv,'%')});
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     %CURVATURE FROM DS 1
%     numNeighbours = 9;
%     scatter_width = 3;%15;
%     clim_set = [0,0.01];
%     gridStep = round(dist*100)/100*2; 
%     ptCloud = pointCloud(PClast); 
%     ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);
%     PClastds = ptCloudA.Location;
%     logicNaN = isnan(PClastds(:,1));
%     idx1 = find(~logicNaN);
%     [normals,curvature] = findPointNormals(PClastds(idx1,:),numNeighbours,[0,0,100],true);
%     temp = find(curvature > prctile(curvature,ptile_curv)); 
%     mask_curv = false(size(PClastds,1),1); 
%     mask_curv(idx1(temp)) = true; 
%     clim_set = [0,prctile(curvature,ptile_curv)];
%     idx3 = find(PClast_mask_in2); 
%     
%     PClastds_mask_in = fit_xyplane_below_threshold_pc(PClastds,guess_threshold,'margin',2*margin_around_plane);
%     
%     %Look at the distance between neighborhood points and detect
%     %outliers according to a multiple (threshold) of the variance.
%     [~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(PClast(PClastds_mask_in&~mask_curv,:)),...
%         'NumNeighbors',10,...
%         'Threshold',3);
%     temp = false(length(inlierIndices),1);
%     temp(inlierIndices) = true;
%     idx_in = find(PClastds_mask_in&~mask_curv);
%     PClastds_mask_in = false(size(PClastds_mask_in)); 
%     PClastds_mask_in(idx_in(temp)) = true; 
%     
%     
%     ca_ds1 = findboundary(PClastds(PClastds_mask_in&~mask_curv,:));
%     
%     % fit ellipse from last contact contour
%     [~,centerEl] = myfitellipse(ca_ds1(:,1),ca_ds1(:,2)); 
%     maskfront = ca_ds1(:,2)<centerEl(2); 
%     [cafit_ds1,~] = myfitellipse(ca_ds1(:,1),ca_ds1(:,2)); 
%     cafit_ds1 = [cafit_ds1,repelem(mean(ca(:,3)),size(cafit_ds1,1),1)]; 
%     
%     subplot(2,4,3); hold on; title(sprintf('Curvature from downsampled pt\ngridstep = %1.3fmm',gridStep));
%     scatter3(PClastds(idx1,1),PClastds(idx1,2),PClastds(idx1,3),scatter_width,curvature,'filled'); 
% %     scatter3(PClastds(idx1(temp),1),PClastds(idx1(temp),2),PClastds(idx1(temp),3),scatter_width,'r','filled'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     subplot(2,4,3+4); hold on; title(sprintf('Segmentation'));
%     scatter3(PClastds(:,1),PClastds(:,2),PClastds(:,3),scatter_width,curvature,'filled'); 
%     hp1=scatter3(PClastds(PClastds_mask_in&~mask_curv,1),PClastds(PClastds_mask_in&~mask_curv,2),PClastds(PClastds_mask_in&~mask_curv,3),scatter_width,'b','filled'); 
%     hp2=scatter3(PClastds(mask_curv,1),PClastds(mask_curv,2),PClastds(mask_curv,3),scatter_width,'r','filled'); 
%     plot3(ca(:,1),ca(:,2),ca(:,3),'-','LineWidth',1,'Color',[1 1 1].*0.5); 
%     plot3(ca_ds1(:,1),ca_ds1(:,2),ca_ds1(:,3),'k-','LineWidth',3);
% %     plot3(cafit_ds1(:,1),cafit_ds1(:,2),cafit_ds1(:,3),'k-','LineWidth',2,'Color','y');  
%     legend([hp1,hp2],{'in',sprintf('curv>%d%s',ptile_curv,'%')});
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     %CURVATURE FROM DS 2
%     numNeighbours = 9;
%     scatter_width = 6;%15;
%     clim_set = [0,0.01];
%     gridStep = round(dist*100)/100*3; 
%     ptCloud = pointCloud(PClast); 
%     ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);
%     PClastds = ptCloudA.Location;
%     logicNaN = isnan(PClastds(:,1));
%     idx1 = find(~logicNaN);
%     mask_below_threshold = PClast(idx1,3) < -1; 
%     [normals,curvature] = findPointNormals(PClastds(idx1,:),numNeighbours,[0,0,100],true);
%     idx2 = find(curvature > prctile(curvature,ptile_curv)); 
%     mask_curv = false(size(PClastds,1),1); 
%     mask_curv(idx1(idx2)) = true; 
%     clim_set = [0,prctile(curvature,ptile_curv)];
%     idx3 = find(PClast_mask_in2); 
%     
%     PClastds_mask_in = fit_xyplane_below_threshold_pc(PClastds,guess_threshold,'margin',2*margin_around_plane);
%     
%     %Look at the distance between neighborhood points and detect
%     %outliers according to a multiple (threshold) of the variance.
%     [~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(PClast(PClastds_mask_in&~mask_curv,:)),...
%         'NumNeighbors',7,...
%         'Threshold',3);
%     temp = false(length(inlierIndices),1);
%     temp(inlierIndices) = true;
%     idx_in = find(PClastds_mask_in&~mask_curv);
%     PClastds_mask_in = false(size(PClastds_mask_in)); 
%     PClastds_mask_in(idx_in(temp)) = true; 
%     
%     ca_ds2 = findboundary(PClastds(PClastds_mask_in&~mask_curv,:));
%     
%     % fit ellipse from last contact contour
%     [~,centerEl] = myfitellipse(ca_ds2(:,1),ca_ds2(:,2)); 
%     maskfront = ca_ds2(:,2)<centerEl(2); 
%     [cafit_ds2,~] = myfitellipse(ca_ds2(maskfront,1),ca_ds2(maskfront,2)); 
%     cafit_ds2 = [cafit_ds2,repelem(mean(ca(:,3)),size(cafit_ds2,1),1)]; 
%     
%     subplot(2,4,4); hold on; title(sprintf('Curvature from downsampled pt\ngridstep = %1.3fmm',gridStep));
%     scatter3(PClastds(idx1,1),PClastds(idx1,2),PClastds(idx1,3),scatter_width,curvature,'filled'); 
% %     scatter3(PClastds(idx1(idx2),1),PClastds(idx1(idx2),2),PClastds(idx1(idx2),3),scatter_width,'r','filled');
%     
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     subplot(2,4,4+4); hold on; title(sprintf('Segmentation'));
%     scatter3(PClastds(:,1),PClastds(:,2),PClastds(:,3),scatter_width,curvature,'filled'); 
%     hp1=scatter3(PClastds(PClastds_mask_in&~mask_curv,1),PClastds(PClastds_mask_in&~mask_curv,2),PClastds(PClastds_mask_in&~mask_curv,3),scatter_width,'b','filled'); 
%     hp2=scatter3(PClastds(mask_curv,1),PClastds(mask_curv,2),PClastds(mask_curv,3),scatter_width,'r','filled'); 
%     plot3(ca(:,1),ca(:,2),ca(:,3),'-','LineWidth',1,'Color',[1 1 1].*0.5); 
%     plot3(ca_ds1(:,1),ca_ds1(:,2),ca_ds1(:,3),'k-','LineWidth',3);
% %     plot3(cafit_ds1(:,1),cafit_ds1(:,2),cafit_ds1(:,3),'-','LineWidth',2,'Color','y');  
%     legend([hp1,hp2],{'in',sprintf('curv>%d%s',ptile_curv,"%")});
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     subplot(2,4,4); hold on; 
%     scatter3(PClast(idx1,1),PClast(idx1,2),PClast(idx1,3),scatter_width,curvature,'filled'); 
%     scatter3(PClast(PClast_mask_in3,1),PClast(PClast_mask_in3,2),PClast(PClast_mask_in3,3),scatter_width,'b','filled'); 
%     scatter3(PClast(mask_curv,1),PClast(mask_curv,2),PClast(mask_curv,3),scatter_width,'r','filled'); 
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,0]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     [xyFitCurv] = myfitellipse(PClastds(idx1(idx2),1),PClastds(idx1(idx2),2)); 
%     
%     scatter3(PClastds(~logicNaN,1),PClastds(~logicNaN,2),PClastds(~logicNaN,3),scatter_width,curvature,'filled'); 
%     
%     idx_IN_curv = inpolygon(PClast(:,1),PClast(:,2),xyFitCurv(:,1),xyFitCurv(:,2));
%     
% %     notNaN = ~logicNaN; 
% %     scatter3(PClast(:,1),PClast(:,2),PClast(:,3),1,'filled'); 
%     scatter3(PClastds(~logicNaN,1),PClastds(~logicNaN,2),PClastds(~logicNaN,3),scatter_width,curvature,'filled'); 
% %     scatter3(PClast(temp3,1),PClast(temp3,2),PClast(temp3,3),1,'r');
%     scatter3(PClast(PClast_mask_in3,1),PClast(PClast_mask_in3,2),PClast(PClast_mask_in3,3),2,'b','filled');
%     scatter3(PClast(idx_IN_curv,1),PClast(idx_IN_curv,2),PClast(idx_IN_curv,3),4,'r','filled');
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'clim',clim_set)
%     colormap summer
%     colorbar
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     hold off;
%     subplot(1,3,3); boxplot(curvature); yline(prctile(curvature,75))





%     figure; hold on; 
%     scatter3(PCnow(:,1),PCnow(:,2),PCnow(:,3),1,[1 1 1].*0.5);
%     plot3(PCnow(PCnow_in_border_id,1),PCnow(PCnow_in_border_id,2),PCnow(PCnow_in_border_id,3),'r');
%     plot3(CAnow(:,1),CAnow(:,2),CAnow(:,3),'b','LineWidth',2); 
%     
%     xlim(gca,[-15,15]); 
%     ylim(gca,[-10,25]); 
%     zlim(gca,[-2,15]); 
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     hold off; 




%         fig = figure('visible','off'); hold on; ax = gca; 
%         title(ax_set,sprintf('Curvature\n from original pt'));
%         scatter3(ax_set,PCnow(:,1),PCnow(:,2),PCnow(:,3),2,curvature_now,'filled'); 
%         xlim(ax_set,[-15,15]);
%         ylim(ax_set,[-10,25]);
%         zlim(ax_set,[-2,15]);
%         set(ax_set,'clim',[0 threshold_absolute_last])
%         colormap(ax_set,'summer'); 
% %         colorbar
%         set(ax_set,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
        
%         format = 'png'; 
%         ii = 1; 
%         fileNameNow = sprintf('curv_%03d.%s',ii,format); 
%         while exist(fileNameNow,'file') 
%             ii = ii + 1; 
%             fileNameNow = sprintf('curv_%03d.%s',ii,format); 
%         end
%         fileNameNow
%         exportgraphics(ax,fileNameNow,...
%             'BackgroundColor','none',...
%             'ContentType','vector');



%             % first elements boundary 
%             ca_1 = findboundary(PCnow(mask_in,:),'connex_factor',0.5);            
%             % get indices of first border
%             [~,idx_ca_1] = ismember(ca_1,PCnow,'rows');
%             mask_ca_1 = false(size(mask_in)); 
%             mask_ca_1(idx_ca_1) = true; 
%             % first elements boundary 
%             ca_2 = findboundary(PCnow(mask_in&~mask_ca_1,:),'connex_factor',0.5);
            % concatenate two successive boundary to get a more meaningful
            % ellipse fitting. 
%             CAnow = cat(1,ca_1,ca_2); 
%             CAnow = ca_1; 
            %Not enough points in contact to differenciate the contact
            %using curvature as proxy (information)
            % first elements boundary
             % fit ellipse from last contact contour : 
        % We assume perfect symmetry and modify the most proximal part of
        % the ellipse to be the exact mirror of the distal part (part where 
        % should be the most confident about measurement). 
%         [~,centerEl] = myfitellipse(ca(:,1),ca(:,2)); 
% %         maskfront = ca(:,2)<centerEl(2); 
% %         ca_front = ca(maskfront,1:2); 
% %         %mirror coordinate  
% %         ca_bottom = [ca(maskfront,1),2*centerEl(2)-ca(maskfront,2)]; 
% %         %concatenation
% %         ca_modif = [ca_front;ca_bottom]; 
% %         %fit
% %         [cafit,~] = myfitellipse(ca_modif(:,1),ca_modif(:,2)); 
%         [cafit,~] = myfitellipse(ca(:,1),ca(:,2));
%         CAnow = [cafit,repelem(mean(ca(:,3)),size(cafit,1),1)]; 
        
%         if high_bound_threshold>mean(PCnow(PCnow_in_id,3)) 
%             temp = PCnow; 
%             temp(~PCnow_in_id,:) = NaN; 
%             CAnow = findboundary(temp);
%         else
%             PCnow_in_id = false(size(PCnow,1),1);  
%         end