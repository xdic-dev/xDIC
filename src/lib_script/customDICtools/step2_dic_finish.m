
function step2_dic_finish(basePath, varargin)

% Define input parser and validation rules
p = inputParser;
p.addRequired('basePath', @isText);
p.addParameter('cam_1',[], @isnumeric);
p.addParameter('cam_2',[], @isnumeric);

p.parse(basePath, varargin{:});
cam_1 = p.Results.cam_1;
cam_2 = p.Results.cam_2;

%% data loading
ncorr1 = load(fullfile(basePath, "ncorr"+sprintf("%d",cam_1)+".mat"));
ncorr2 = load(fullfile(basePath, "ncorr"+sprintf("%d",cam_2)+".mat"));
ncorr1_2 = load(fullfile(basePath, "ncorr"+sprintf("%d",[cam_1,cam_2])+".mat"));

data1 = ncorr1.data_dic_save;
data1_2 = ncorr1_2.data_dic_save;
data2 = ncorr2.data_dic_save;

% disp('loading Done');

%% filtering
% data1 = test_filtering_field(data1); 
% data1_2 = test_filtering_field(data1_2); 
% data2 = test_filtering_field(data2); 
% disp('cleaning Done');

%% spatio-temporal low pass filter to displacement field on subset with low 
% correlation coefficent.
fig = figure('Visible','off'); h_ax = subplot_ax(3,1); 
data1 = replacebadcorr(data1,h_ax(1));
data1_2 = replacebadcorr(data1_2,h_ax(2));
data2 = replacebadcorr(data2,h_ax(3));
saveas(fig,fullfile(basePath,sprintf('disp_filt_%d-%d.png',cam_1,cam_2))); 

%% fill void
% window = 3; 
% data1 = fill_void(data1, window); 
% data2 = fill_void(data2, window);
% data1_2 = fill_void(data1_2, window);

%% start DIC data structure
DIC_results= struct;

DIC_results.nCamRef = cam_1;
DIC_results.nCamDef = cam_2;
DIC_results.nImages = size(data1.displacements,2);

DIC_results.ROImask = ncorr1.reference_save.roi.mask;

DIC_results.ncorrInfo = data1.dispinfo;
DIC_results.ncorrInfo.cutoff_corrcoef = [
    DIC_results.ncorrInfo.cutoff_corrcoef;
    data1_2.dispinfo.cutoff_corrcoef;
    data2.dispinfo.cutoff_corrcoef;
    ];

% disp('init Done')

%% compute points / corCoeffVec / Faces / FaceColors

% Extract results
Disp1=data1.displacements;
DispInfo1=data1.dispinfo;

Disp1_2=data1_2.displacements;
DispInfo1_2=data1_2.dispinfo;

Disp2=data2.displacements;
DispInfo2=data2.dispinfo;

nCur=size(Disp1,2) + 1 + size(Disp2,2);

Factor=DispInfo1.spacing+1;

% Define output structure
ROI_DIC=cell(nCur,1);
CorCoeff=cell(nCur,1); 
CorCoeffVec=cell(nCur,1); 
Points=cell(nCur,1);
Uvec=cell(nCur,1); 
Vvec=cell(nCur,1);


[YrefROIVec1,XrefROIVec1] = find(Disp1(1).roi_dic.mask);
PtempRef1=[XrefROIVec1,YrefROIVec1];
PtempRef1=(PtempRef1-1)*Factor+1; % switch from spacing to pixels
Pref=PtempRef1;

IMref = ncorr1.reference_save.gs; 

% disp('data ready')

%% First loop on data1
outside_it = 0;
for ii=1:size(Disp1,2)

    outside_it = outside_it+1;
    
    ROI_DIC{ii}=Disp1(ii).roi_dic.mask;
    
    CorCoeff{ii}=Disp1(ii).plot_corrcoef_dic;
    Uref=Disp1(ii).plot_u_ref_formatted;
    Vref=Disp1(ii).plot_v_ref_formatted;
       
    CorCoeffVec{ii}=CorCoeff{ii}(ROI_DIC{1});
    CorCoeffVec{ii}(CorCoeffVec{ii}==0)=NaN;
    
    % displacements from ref to cur
    UrefROIVec=Uref(ROI_DIC{1});
    UrefROIVec(UrefROIVec==0)=NaN;
    VrefROIVec=Vref(ROI_DIC{1});
    VrefROIVec(VrefROIVec==0)=NaN;
    
    Uvec{ii}=UrefROIVec;
    Vvec{ii}=UrefROIVec;
    % current points
    Points{ii}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];
    
    % save face colors for further 3D analysis
    if ii==1
        % pixel colors
        IMrefSmall=IMref(1:Factor:end,1:Factor:end);
        IMrefSmallMasked=IMrefSmall;
        IMrefSmallMasked(~ROI_DIC{1})=[];
        ColorRef=IMrefSmallMasked(:);
    end
end

% disp('First loop done')

%% Match between ref of the 2 sets

outside_it=outside_it+1;
ii=outside_it;

ROI_DIC{ii}=Disp1_2(1).roi_dic.mask;

CorCoeff{ii}=Disp1_2(1).plot_corrcoef_dic;
Uref=Disp1_2(1).plot_u_ref_formatted;
Vref=Disp1_2(1).plot_v_ref_formatted;
   
U1_2 = Uref;
V1_2 = Vref;

CorCoeffVec{ii}=CorCoeff{ii}(ROI_DIC{1});
CorCoeffVec{ii}(CorCoeffVec{ii}==0)=NaN;

% displacements from ref to cur
UrefROIVec=Uref(ROI_DIC{1});
UrefROIVec(UrefROIVec==0)=NaN;
VrefROIVec=Vref(ROI_DIC{1});
VrefROIVec(VrefROIVec==0)=NaN;

Uvec{ii}=UrefROIVec;
Vvec{ii}=UrefROIVec;

% current points
Points{ii}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];

% disp('Matching 1-2 done')

%% compute point mapping
U_mapped = U1_2/Factor;
V_mapped = V1_2/Factor;

[U_mapped_row, U_mapped_col, U_mapped_value] = find(U_mapped);
[V_mapped_row, V_mapped_col, V_mapped_value] = find(V_mapped);

mapped_indices = sub2ind(size(U_mapped), U_mapped_row, U_mapped_col);

U_down = floor(U_mapped_value);
U_up = ceil(U_mapped_value);
V_down = floor(V_mapped_value);
V_up = ceil(V_mapped_value);

U_low = U_mapped_col + U_down;
U_high = U_mapped_col + U_up;
V_low = V_mapped_row + V_down;
V_high = V_mapped_row + V_up;

ind_low_low = sub2ind(size(U_mapped), V_low, U_low);
ind_high_low = sub2ind(size(U_mapped), V_low, U_high);
ind_low_high = sub2ind(size(U_mapped), V_high, U_low);
ind_high_high = sub2ind(size(U_mapped), V_high, U_high);

real_col = U_mapped_col + U_mapped_value;
real_row = V_mapped_row + V_mapped_value;

% compute distances
s(:,1) = 1./sqrt((real_col-U_low).^2 + (real_row-V_low).^2);
s(:,2) = 1./sqrt((real_col-U_high).^2 + (real_row-V_low).^2);
s(:,3) = 1./sqrt((real_col-U_low).^2 + (real_row-V_high).^2);
s(:,4) = 1./sqrt((real_col-U_high).^2 + (real_row-V_high).^2);


%% Second loop on data2
for ii=1:size(Disp2,2)

    outside_it = outside_it+1;
    
    ROI_DIC{outside_it}=Disp2(ii).roi_dic.mask;
    
    CorCoeff{outside_it}=Disp2(ii).plot_corrcoef_dic;
    Uref=Disp2(ii).plot_u_ref_formatted;
    Vref=Disp2(ii).plot_v_ref_formatted;

    all_u = [
        Uref(ind_low_low) Uref(ind_high_low) Uref(ind_low_high) Uref(ind_high_high)
        ];

    all_v = [Vref(ind_low_low) Vref(ind_high_low) Vref(ind_low_high) Vref(ind_high_high)
        ];

    scaled_u = sum(all_u.*s, 2)./sum((all_u~=0).*s, 2);
    scaled_v = sum(all_v.*s, 2)./sum((all_v~=0).*s, 2);
    
    scaled_u(scaled_u==0)=nan;
    scaled_v(scaled_v==0)=nan;

    new_Uref = zeros(size(Uref));
    new_Vref = zeros(size(Vref));
    
    new_Uref(mapped_indices) = scaled_u;
    new_Vref(mapped_indices) = scaled_v;

    % correct the displacement
    Uref = new_Uref + U1_2;
    Vref = new_Vref + V1_2;
       
    CorCoeffVec{outside_it}=CorCoeff{outside_it}(ROI_DIC{1});
    CorCoeffVec{outside_it}(CorCoeffVec{outside_it}==0)=NaN;
    
    % displacements from ref to cur
    UrefROIVec=Uref(ROI_DIC{1});
    UrefROIVec(UrefROIVec==0)=NaN;
    VrefROIVec=Vref(ROI_DIC{1});
    VrefROIVec(VrefROIVec==0)=NaN;
    
    Uvec{outside_it}=UrefROIVec;
    Vvec{outside_it}=UrefROIVec;
    % current points
    Points{outside_it}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];    
end

% disp('Second loop done')

%% Finish step 2
% Create triangulation
DT = delaunayTriangulation(Pref);
F = DT.ConnectivityList;
V = Pref;

% remove irregular triangles
EdgeLengths = patchEdgeLengths(F,V);
EdgeLengths = [EdgeLengths(1:3:length(EdgeLengths)) EdgeLengths(2:3:length(EdgeLengths)) EdgeLengths(3:3:length(EdgeLengths))];
EdgeLengthsMax=max(EdgeLengths,[],2);

F(EdgeLengthsMax>1.1*sqrt(2)*Factor,:)=[];

% flip direction to have normals pointing out
F=F(:,[1 3 2]);

% face colors (average node colors)
CF=mean(ColorRef(F),2);

% disp('Triangulation Done');

%% Add data to output structure & save
DIC_results.Points=Points;
DIC_results.CorCoeffVec=CorCoeffVec;
DIC_results.Faces=F;
DIC_results.FaceColors=CF;


output_filename = fullfile( basePath, ...
    sprintf("myDIC2DpairResults_C_%d_C_%d.mat", cam_1, cam_2));

DIC2DpairResults = DIC_results;
save(output_filename, "DIC2DpairResults");

% disp('File Saved - Finished')

end


% function output = isText(myVar)
% % isText checks if the input variable is a text type.
% % The function returns true if the variable is of string, character array, or cell array of strings type; otherwise, it returns false.
% %
% % Input:
% %   - myVar: The variable to be checked.
% %
% % Output:
% %   - output: A logical value indicating if the variable is a text type.
% 
%     output = isstring(myVar) || ischar(myVar) || iscellstr(myVar);
% end
% 
% function data = fill_void(data, window)
% Factor = data.dispinfo.spacing + 1; 
% Nframe = length(data.displacements); 
% siz = size(data.displacements(1).plot_u_ref_formatted); 
% u = zeros(siz(1),siz(2),size(data.displacements,2));
% v = zeros(siz(1),siz(2),size(data.displacements,2));
% mask = imresize(data.displacements(1).roi_dic.mask,1/Factor);
% [y,x] = find(mask== 1); 
% idx2clean = cell(Nframe,1); 
% points2clean = zeros(length(y),2,Nframe); 
% points2clean_cell = cell(length(y),1);
% 
% for ii = 1:Nframe
%     u_ii = data.displacements(ii).plot_u_ref_formatted;
%     v_ii = data.displacements(ii).plot_v_ref_formatted;
%     u(:,:,ii) = u_ii; 
%     v(:,:,ii) = v_ii; 
%     idx2clean{ii} = find(u_ii(mask) == 0);
%     points2clean(1:length(idx2clean{ii}),:,ii) = [x(idx2clean{ii}),y(idx2clean{ii})]; 
%     points2clean_cell{ii} = [x(idx2clean{ii}),y(idx2clean{ii})]; 
% end
% ufilled = fill_displ_void(u, points2clean_cell, window);
% vfilled = fill_displ_void(v, points2clean_cell, window);
% 
% for ii = 1:Nframe
%     data.displacements(ii).plot_u_ref_formatted = ufilled(:,:,ii); 
%     data.displacements(ii).plot_v_ref_formatted = vfilled(:,:,ii); 
% end
% end
% 
% function dispfilled = fill_displ_void(disp, points2clean_cell, window)
%     dispfilled = disp; 
%     for ii = 1:size(disp,3)
%         pii = points2clean_cell{ii};
%         L = length(pii); 
%         for hh = 1:L
%             disparound = disp(pii(hh,2)+[-window:window],pii(hh,1)+[-window:window],ii);
%             disp_mean = nanmean(disparound(disparound~=0));
%             disp(size(disp_mean))
%             if ~isnan(disp_mean)
%                 dispfilled(pii(hh,2),pii(hh,1),ii) = disp_mean;
%             end
%         end
%     end
% end
% function data = cleandisplacementfield(data,im,plotting)
%     
%     movemean = 15; 
%     
%     Nframe = size(data.displacements,2);  
%     siz = size(data.displacements(1).plot_u_ref_formatted); 
%     u = zeros(siz(1),siz(2),size(data.displacements,2));
%     v = zeros(siz(1),siz(2),size(data.displacements,2));
%     for i = 1:Nframe
%         u(:,:,i) = data.displacements(i).plot_u_ref_formatted;
%         v(:,:,i) = data.displacements(i).plot_v_ref_formatted;
%     end
% 
%     %outlier filling method on the whole sequence is performed for the sum of 
%     % displacement field to find outlier and interpolate in each displacement 
%     % field.  
%     [ufilt,vfilt,outlier] = myoutlierfilling(u,v,movemean);
%     if plotting
%         if ~isempty(outlier)
%             figure; hold on;
%             for i = 1:size(outlier,1)
%                 plot(1:Nframe,squeeze(u(outlier(i,1),outlier(i,2),:)+v(outlier(i,1),outlier(i,2),:)),'-','Color',[1 1 1].*0.8);
%                 plot(1:Nframe,squeeze(ufilt(outlier(i,1),outlier(i,2),:)+vfilt(outlier(i,1),outlier(i,2),:)),'k--');
%                 h1 = plot(outlier(i,3),squeeze(u(outlier(i,1),outlier(i,2),outlier(i,3))+...
%                     v(outlier(i,1),outlier(i,2),outlier(i,3))),'rx');
%                 h2 = plot(outlier(i,3),squeeze(ufilt(outlier(i,1),outlier(i,2),outlier(i,3))+...
%                     vfilt(outlier(i,1),outlier(i,2),outlier(i,3))),'bx');
%             end
%             legend([h1,h2],'original','replaced'); 
%             set(gca,'Color','None'); 
%             
%             %
%             points = zeros(size(outlier,1),2,Nframe); 
%             pointsfilt = zeros(size(outlier,1),2,Nframe); 
%             factor = 0.1; 
%             points(:,1,1) = outlier(:,2); 
%             points(:,2,1) = outlier(:,1); 
%             pointsfilt(:,1,1) = outlier(:,2); 
%             pointsfilt(:,2,1) = outlier(:,1); 
%             imr = zeros(siz(1),siz(2),Nframe,'uint8'); 
%             imr(:,:,1) = imresize(im(:,:,1),factor,'nearest');
%             for i = 1:Nframe-1
%                 points(:,1,i+1) = (points(:,1,1)+factor*diag(squeeze(u(outlier(:,1)', outlier(:,2)',i)))); 
%                 points(:,2,i+1) = (points(:,2,1)+factor*diag(squeeze(v(outlier(:,1)', outlier(:,2)',i)))); 
%                 pointsfilt(:,1,i+1) = (pointsfilt(:,1,1)+factor*diag(squeeze(ufilt(outlier(:,1)', outlier(:,2)',i)))); 
%                 pointsfilt(:,2,i+1) = (pointsfilt(:,2,1)+factor*diag(squeeze(vfilt(outlier(:,1)', outlier(:,2)',i)))); 
%                 imr(:,:,i+1) = imresize(im(:,:,i+1),factor,'nearest');
%             end
% 
%             player(imr,'points',points,'morepoints',pointsfilt);
%         end
%     end
%     for i = 1:Nframe
%         data.displacements(i).plot_u_ref_formatted = ufilt(:,:,i);
%         data.displacements(i).plot_v_ref_formatted = vfilt(:,:,i);
%     end
%     
% end
% function [ufilt,vfilt,outlier] = myoutlierfilling(u,v,movmeanvalue)
%     siz = size(u);
%     u = permute(u, [3 1 2]);
%     v = permute(v, [3 1 2]);
%     ufilt = zeros(siz(3),siz(1),siz(2)); vfilt = zeros(siz(3),siz(1),siz(2));
%     x=[]; y=[]; img=[];
%     for i = 1:siz(2) % loop on the original column of the input matrix.
%         %find outliers from the global displacement field
%         idx=isoutlier(squeeze(u(:,:,i)+v(:,:,i)),'movmean',movmeanvalue);
%         xi = []; 
%         [imgi,xi] = find(idx==1);
%         if ~isempty(xi)
%             x = [x;xi]; 
%             y = [y;repelem(i,length(xi))'];
%             img = [img;imgi];
%         end
%         %fill outliers
%         ufilt(:,:,i) = filloutliers(u(:,:,i),'linear','OutlierLocations', idx);
%         vfilt(:,:,i) = filloutliers(v(:,:,i),'linear','OutlierLocations', idx);
%     end
%     ufilt = permute(ufilt, [2 3 1]);
%     vfilt = permute(vfilt, [2 3 1]);
%     outlier = [x,y,img];
% end 

 function [handles] = replacebadcorr(handles,h_ax_ii)
    
    %Replacement of all badly correlated data by filtered data
    level_corr_coef = 1; 

    Nframe = length(handles.displacements);
    siz = size(handles.displacements(end).plot_u_ref_formatted);
    if Nframe == 2 % matching step
        Nframe = 1; 
    end
    
    %roi 
    M = imresize(handles.displacements(end).plot_u_ref_formatted,[siz(1),siz(2)]); 
    
    %initialization
    z = zeros(siz(1),siz(2),Nframe); 
    u = z; v = z; corr_coef = z;
    
    for ii = 1:Nframe
        u_ii = handles.displacements(ii).plot_u_ref_formatted; 
        v_ii = handles.displacements(ii).plot_v_ref_formatted; 
        corr_coef_ii=handles.displacements(ii).plot_corrcoef_dic; %metric
        
%         [u_ii,~] = fillmissing2(u_ii,"cubic",MissingLocations=M); 
%         [v_ii,~] = fillmissing2(v_ii,"cubic",MissingLocations=M); 
%         [corr_coef_ii,~] = fillmissing2(corr_coef_ii,"cubic",MissingLocations=M); 
        
        u(:,:,ii) = u_ii; 
        v(:,:,ii) = v_ii; 
        corr_coef(:,:,ii) = corr_coef_ii;
    end
    mask_corr_coef = logical(corr_coef>level_corr_coef); 
    umodif1 = u; vmodif1 = v; 
    umodif2 = u; vmodif2 = v; 
    
    %initial filtered maps with badly correlated points : 
    %replicates is used in order to minimize the border effect
    sigma = [1 1 2]; 
    uf = imgaussfilt3(u, sigma,'FilterDomain','spatial','padding','replicate');
    vf = imgaussfilt3(v, sigma,'FilterDomain','spatial','padding','replicate');

    % remove border effect of the spatial filtering on subset not very well
    % correlated.
    mask_in = u(:,:,end)~=0; 
    [mask_out] = shrink_mask(mask_in,0.9);
    mask_out=repmat(mask_out,[1 1 Nframe]);

    umodif1(mask_corr_coef&mask_out)=uf(mask_corr_coef&mask_out);
    vmodif1(mask_corr_coef&mask_out)=vf(mask_corr_coef&mask_out);
    
    %second filtered maps with badly correlated points replaced :
    %this step makes a smoother filtered map around the badly correlated
    %points, or said differently, suppress more of the supposedly bad 
    %influence of this points.
    uf = imgaussfilt3(umodif1, sigma,'FilterDomain','spatial','padding','replicate');
    vf = imgaussfilt3(vmodif1, sigma,'FilterDomain','spatial','padding','replicate');
    
    umodif2(mask_corr_coef&mask_out)=uf(mask_corr_coef&mask_out);
    vmodif2(mask_corr_coef&mask_out)=vf(mask_corr_coef&mask_out);
    
    sub_repl = sum(mask_corr_coef&mask_out,'all','omitnan'); 
    fprintf('nbr replaced subset : %d\n',sub_repl); 
    
    set(gcf,'CurrentAxes',h_ax_ii);
    imshow(cat(1,cat(2,u(:,:,end),uf(:,:,end),umodif2(:,:,end)),...
                 cat(2,v(:,:,end),vf(:,:,end),vmodif2(:,:,end))),[]); 
    colormap(h_ax_ii,'jet'); 
    title(h_ax_ii,sprintf('replaced : %d',sub_repl));
             
    for ii = 1:Nframe
        handles.displacements(ii).plot_u_ref_formatted = umodif2(:,:,ii); 
        handles.displacements(ii).plot_v_ref_formatted = vmodif2(:,:,ii); 
    end
 end

 
 
 
 
 