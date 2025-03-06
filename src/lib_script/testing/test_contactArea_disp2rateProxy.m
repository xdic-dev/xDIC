
r = DIC3DPPresults; 
%%
xyz=r.FaceCentroids;
[~,DispZmat_rate] = data_dic2matrix(r,'DispZ'); 
%%
figure; 
colormap jet 
colorbar 
set(gca,'clim',[-0.15 0.15])
set(gca,'Color','none'); 
%%

%Outliers detection on individual maps and 
%Combine individual outliers id. 
mask_outliers_indiv_id = cell(1,3);
id = [1,36+10,length(r.FaceCentroids)];
for iframeMap = 1:length(id)
    ptFC_cur = xyz{id(iframeMap)};

    isNotValid = isnan(ptFC_cur(:,1));

    %Look at the distance between neighborhood points and detect
    %outliers according to a multiple (threshold) of the variance.
    [~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(ptFC_cur),...
        'NumNeighbors',15,...
        'Threshold',3);
    temp = false(length(ptFC_cur),1);
    temp(outlierIndices) = true;
    mask_outliers_indiv_id{iframeMap} = temp|isNotValid;
end
mask_comb_current = false(length(r.FaceCentroids{1}),1);
for iframeMap = 1:length(id)
    mask_comb_current = mask_comb_current|mask_outliers_indiv_id{iframeMap}; 
end

%%
[nPoint, nFrame] = size(DispZmat_rate); 
contact_start = 27; 
id = [contact_start-2:contact_start+60]; 
theshold_percent = 75; 
xyz_last = xyz{77};
min_point_inside = 3; 
margin_around_plane = 0.1; 
color_scale = jet(length(id)); %1-linspace(0.1,0.9,length(id));    

newfig('test_dispZrate_proxy2'); 
nRow = 2; 
nCol = 4;%length(id); 
merge_element = {}; %{(nRow-1)*nCol+1:nRow*nCol};
ax = subplot_ax(nRow,nCol,'merge',merge_element); hold(ax,'on'); 
disp_counter = 0; 
for ii = 1:length(id)
    color_scale_ii = color_scale(ii,:); 
    iframe = id(ii); 
    xyz_now = xyz{iframe}(:,:);
    xyz_now(mask_comb_current,:) = NaN;
    DispZmat_rate_now =abs(DispZmat_rate(:,iframe)); 
    DispZmat_rate_now(mask_comb_current,:) = NaN;
    
    mask_sub = -1>xyz_now(:,3);  
    idx_sub = find(mask_sub); 
    max_disp_current = prctile(DispZmat_rate_now(mask_sub,:),95,1); 
    actual_threshold = prctile(DispZmat_rate_now(mask_sub,:),theshold_percent,1); 
    actual_threshold = 0.01; 
    
    mask_in = actual_threshold>DispZmat_rate_now(:,:); 
    temp = find(mask_in);
    idx_in = temp; %idx_sub(temp); 
    
    CAnow1 = [NaN,NaN,NaN];
    disp_logic = false; 
    if length(idx_in) > min_point_inside %&& max_disp_current > actual_threshold*5
        [~,inlierIndices,~] = pcdenoise(pointCloud(xyz_now(idx_in,:)),...
            'NumNeighbors',20,...
            'Threshold',3);
        temp = false(length(inlierIndices),1);
        temp(inlierIndices) = true;

        idx_in = find(mask_in);
        mask_in = false(size(mask_in)); 
        mask_in(idx_in(temp)) = true; 
        CAnow1 = findboundary(xyz_now(mask_in & mask_sub,:));
%         CAnow1 = [CAnow1,repelem(-2,length(ca)),1];
        disp_logic = true; 
    end
    
%     ax_idx = ii; 
%     scatter3(ax(ax_idx),xyz_now(:,1),xyz_now(:,2),xyz_now(:,3),1,[1 1 1].*0.5); 
%     scatter3(ax(ax_idx),xyz_now(idx_in,1),xyz_now(idx_in,2),xyz_now(idx_in,3),5,'r'); 
%     plot3(ax(ax_idx),CAnow(:,1),CAnow(:,2),CAnow(:,3)-0.1,'k','LineWidth',3); 
%     xlim(ax(ax_idx),[-10,10]);
%     ylim(ax(ax_idx),[-5,15]);
%     zlim(ax(ax_idx),[-2,15]);
%     set(ax(ax_idx),'clim',[-0.15 0.15])
%     colormap(ax(ax_idx),'summer'); 
%     set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     ax_idx = ii+nCol; 
%     scatter3(ax(ax_idx),xyz_now(:,1),xyz_now(:,2),xyz_now(:,3),1,DispZmat_rate_now); 
%     plot3(ax(ax_idx),CAnow(:,1),CAnow(:,2),CAnow(:,3)-0.1,'k','LineWidth',3); 
%     xlim(ax(ax_idx),[-10,10]);
%     ylim(ax(ax_idx),[-5,15]);
%     zlim(ax(ax_idx),[-2,15]);
%     set(ax(ax_idx),'clim',[-0.15 0.15])
%     colormap(ax(ax_idx),'jet'); 
%     set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     if disp_logic && disp_counter < 8
        ax_idx = 1; 
        plot(ax(ax_idx),CAnow1(:,1),CAnow1(:,2),'-','Color',color_scale_ii); 
        xlim(ax(ax_idx),[-10,10]);ylim(ax(ax_idx),[-5,15]);set(ax(ax_idx), 'YDir','reverse'); 
        title(ax(ax_idx),'Proxy z-displacement rate'); 
        [val_min, idx_min] = min(CAnow1(:,1)); 
        text(ax(ax_idx),val_min,CAnow1(idx_min,2),sprintf('%1.2f',nf_all(iframe,1)),...
            'HorizontalAlignment','center','FontSize',12); 
        set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None');
        if ii == length(id) 
            ax_idx = 2; 
            plot(ax(ax_idx),CAnow1(:,1), CAnow1(:,2),'-','Color','k','LineWidth',2); 
            scatter(ax(ax_idx),xyz_now(:,1),xyz_now(:,2),1,'k','filled'); 
            set(ax(ax_idx),'Color','None','DataAspectRatio',[1 1 1],'YDir','reverse','clipping','off')
            xlim(ax(ax_idx),[-10,10]);ylim(ax(ax_idx),[-5,15]);
        end
        
%     end
    
%         ax_idx = 2; 
%         plot(ax(ax_idx),CAnow(:,1),CAnow(:,2),'-','Color',color_scale_ii); 
%         xlim(ax(ax_idx),[-10,10]);ylim(ax(ax_idx),[-5,15]);set(ax(ax_idx), 'YDir','reverse'); 
%         
%         [val_min, idx_min] = min(CAnow(:,1)); 
%         text(ax(ax_idx),val_min,CAnow(idx_min,2),sprintf('%1.2f',nf_all(iframe,1)),...
%             'HorizontalAlignment','center','FontSize',12); 
%         set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None');
%     
    % FROM CURVATURE 
    masknow_below = xyz_now(:,3) < -1; 
    CAnow2 = [NaN,NaN,NaN];
    xyz_now_in_id = false(size(xyz_now,1),1);  

    %CURVATURE FROM ORIGINAL POINT CLOUD
    numNeighbours = 9;
    [~,threshold_absolute_last] = find_curv_threshold(xyz_last,'numNeighbours',numNeighbours,'threshold_percent',50);
    [mask_curv_now,~,curvature_now] = find_curv_threshold(xyz_now,...
        'numNeighbours',numNeighbours,...
        'threshold_percent',30);

    if length(find(masknow_below == 1)) > min_point_inside
        [xyz_now_in_id,fun_planefit] = fit_xyplane_below_threshold_pc(xyz_now,-1,'margin',1*margin_around_plane);
        mask_in = xyz_now_in_id&~mask_curv_now; 

        %find contour
        if sum(mask_in) > 3
            %Look at the distance between neighborhood points and detect
            %outliers according to a multiple (threshold) of the variance.
            [~,inlierIndices,~] = pcdenoise(pointCloud(xyz_last(mask_in,:)),...
                'NumNeighbors',20,...
                'Threshold',3);
            temp = false(length(inlierIndices),1);
            temp(inlierIndices) = true;
            idx_in = find(mask_in);
            mask_in = false(size(mask_in)); 
            mask_in(idx_in(temp)) = true; 
            ca = findboundary(xyz_now(mask_in,:),'connex_factor',0.1); 
        else
            %Not enough points in contact to differenciate the contact with
            %curvature information 
            ca = findboundary(xyz_now(xyz_now_in_id,:),'connex_factor',0.1);
        end
        CAnow2 = [ca];
    end 
    
%     ax_idx = ii + nCol*2;
%     scatter3(ax(ax_idx),xyz_now(:,1),xyz_now(:,2),xyz_now(:,3),1,[1 1 1].*0.5); 
%     scatter3(ax(ax_idx),xyz_now(mask_in,1),xyz_now(mask_in,2),xyz_now(mask_in,3),5,'r'); 
%     plot3(ax(ax_idx),CAnow(:,1),CAnow(:,2),CAnow(:,3),'k','LineWidth',3); 
%     plot3(ax(ax_idx),CAnow(:,1),CAnow(:,2),CAnow(:,3)-0.1,'k','LineWidth',3); 
%     xlim(ax(ax_idx),[-10,10]);
%     ylim(ax(ax_idx),[-5,15]);
%     zlim(ax(ax_idx),[-2,15]);
% %     set(ax(ax_idx),'clim',[-0.15 0.15])
%     colormap(ax(ax_idx),'summer'); 
%     set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     
%     ax_idx = ii + nCol*3; 
%     scatter3(ax(ax_idx),xyz_now(:,1),xyz_now(:,2),xyz_now(:,3),1,mask_curv_now); 
%     plot3(ax(ax_idx),CAnow(:,1),CAnow(:,2),CAnow(:,3)-0.1,'k','LineWidth',3); 
%     xlim(ax(ax_idx),[-10,10]);
%     ylim(ax(ax_idx),[-5,15]);
%     zlim(ax(ax_idx),[-2,15]);
%     set(ax(ax_idx),'clim',[0 threshold_absolute_last])
%     colormap(ax(ax_idx),'jet'); 
%     set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%     if disp_counter >= 8
        ax_idx = 3; 
        plot(ax(ax_idx),CAnow2(:,1),CAnow2(:,2),'-.','Color',color_scale_ii); 
        xlim(ax(ax_idx),[-10,10]);ylim(ax(ax_idx),[-5,15]);set(ax(ax_idx), 'YDir','reverse')
        [val_min, idx_min] = min(CAnow2(:,1)); 
        text(ax(ax_idx),val_min,CAnow2(idx_min,2),sprintf('%1.2f',nf_all(iframe,1)),...
            'HorizontalAlignment','center','FontSize',12); 
        set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None');
        title(ax(ax_idx),'current estim');
        
        if ii == length(id) 
            ax_idx = 4; 
            plot3(ax(ax_idx),CAnow2(:,1), CAnow2(:,2),CAnow2(:,3),'-','Color','k','LineWidth',2); 
            scatter3(ax(ax_idx),xyz_now(:,1),xyz_now(:,2),xyz_now(:,3),1,'k','filled'); 
            set(ax(ax_idx),'View',[0 -90]); 
            set(ax(ax_idx),'Color','None','DataAspectRatio',[1 1 1],'clipping','off')
            xlim(ax(ax_idx),[-10,10]);ylim(ax(ax_idx),[-5,15]);
        end
%     end
    
        ax_idx = 5; 
        disp_counter = disp_counter+1; 
        disp_meth_counter = 10; 
        if disp_counter < disp_meth_counter
            plot(ax(ax_idx),CAnow1(:,1),CAnow1(:,2),'-.','Color',color_scale_ii);
        else
            plot(ax(ax_idx),CAnow2(:,1),CAnow2(:,2),'-.','Color',color_scale_ii);
        end
        xlim(ax(ax_idx),[-10,10]);ylim(ax(ax_idx),[-5,15]);set(ax(ax_idx), 'YDir','reverse')
        [val_min, idx_min] = min(CAnow2(:,1)); 
        text(ax(ax_idx),val_min,CAnow2(idx_min,2),sprintf('%1.2f',nf_all(iframe,1)),...
            'HorizontalAlignment','center','FontSize',12); 
        set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None');
        title(ax(ax_idx),sprintf('meth 1 before %d',disp_meth_counter));
%     ax_idx = ii+2*nCol;
%     boxplot(ax(ax_idx),DispZmat_rate_now(idx_sub)); 
%     axis(ax(ax_idx),[0.8,1.2,0,0.15]); 
end






