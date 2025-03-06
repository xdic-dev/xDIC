function [mask_comb_all,mask_outliers_indiv_id] = point_cloud_outliers_removal(data,id)
    mask_outliers_indiv_id = cell(1,length(id));
    for iframeMap = 1:length(id)
        ptFC_cur = data.FaceCentroids{id(iframeMap)};
        
        isNotValid = isnan(ptFC_cur(:,1));
        
        %Look at the distance between neighborhood points and detect
        %outliers according to a multiple (threshold) of the variance.
        [~,inlierIndices,outlierIndices] = pcdenoise(pointCloud(ptFC_cur),...
            'NumNeighbors',15,...
            'Threshold',3);
        temp = false(length(ptFC_cur),1);
        temp(outlierIndices) = true;
        %     mask_outliers{mm} = temp;
        mask_outliers_indiv_id{iframeMap} = temp|isNotValid;
    end
    mask_comb_all= false(length(data.FaceCentroids{1}),1);
    for iframeMap = 1:length(id)
        mask_comb_all = mask_comb_all|mask_outliers_indiv_id{iframeMap}; 
    end
end