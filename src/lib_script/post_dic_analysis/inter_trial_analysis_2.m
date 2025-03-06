
% Script for post dic analysis 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial selection
SUBJECT = "S10";
PHASE = "loading";
MATERIAL_NBR = 2; 
NFCOND_SET = [5]; 
SPDDXLCOND_SET = [0.04,0.08]; 

% Measurement selection
FACEMEASURE = {'Epc1','Epc2','DispZ','curv'};

% 
trial_target_set = search_trial2target(baseDataPath,SUBJECT,PHASE,FRcond{MATERIAL_NBR},...
    NFCOND_SET,SPDDXLCOND_SET);

%Check that the files exist 
TRIAL_TARGET = []; 
for itrial = 1:length(trial_target_set)
    fileNbrPair = 2; 
    fileversion = "v1"; 
    material = FRcond{MATERIAL_NBR}; 
    file_DIC = fullfile(baseDICPath,SUBJECT,material,sprintf("%03d",trial_target_set(itrial)),PHASE,...
            sprintf('DIC3DPPresults_%dPairs_cum_%s.mat',fileNbrPair,fileversion));
    if exist(file_DIC,'file') == 2
        TRIAL_TARGET = [TRIAL_TARGET,trial_target_set(itrial)];
    end
end

trial_info = search_trialinfo(baseDataPath,SUBJECT,FRcond{MATERIAL_NBR},TRIAL_TARGET);
Ntrial = length(TRIAL_TARGET); 

% Print trials
if length(TRIAL_TARGET)>=2
    fprintf('Analysis of trials [');
    fprintf('%g, ',TRIAL_TARGET(1:end-1)); 
    fprintf('%g]\n',TRIAL_TARGET(end)); 
    fprintf('(withdraw : ');
    fprintf('%g, ',trial_target_set(~ismember(trial_target_set,TRIAL_TARGET)));
    fprintf(')\n');  
else
    fprintf('Analysis of trials [');
    fprintf('%g]\n',TRIAL_TARGET(1)); 
    fprintf('(withdraw : ');
    fprintf('%g, ',trial_target_set(~ismember(trial_target_set,TRIAL_TARGET)));
    fprintf(')\n');  
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Loading trials for subject %s...\n",SUBJECT); 

% Init outputs 
data_dic = cell(Ntrial,1); 
data_robot = cell(Ntrial,1); 

mm = 1;  
for itrial = TRIAL_TARGET
    tic; 
    fprintf("--> loading %d...",itrial);
    fileNbrPair = 2; 
    fileversion = "v1"; 
    material = FRcond{MATERIAL_NBR}; 
    %set filename
    file_DIC = fullfile(baseDICPath,SUBJECT,material,sprintf("%03d",itrial),PHASE,...
        sprintf('DIC3DPPresults_%dPairs_cum_%s.mat',fileNbrPair,fileversion));
    file_robot = fullfile(baseDataPath,'rawdata',SUBJECT,'speckles',material,'robot',...
        sprintf('*_%d.csv',itrial));
    file_protocol = fullfile(fullfile(baseDataPath,'rawdata',SUBJECT,'speckles',material,...
        'protocol',sprintf('*.mat')));
    trialname = SUBJECT+"_"+material+"_speckles_"+sprintf("%03d",itrial);
    % file_cam = fullfile(baseDataPath,'rawdata',subject,"speckles",material,"vid",...
    %     sprintf('%s*_cam_%d*.mp4',trialname,camNbr));
    try
        %dic data
        data_dic_cur = load(file_DIC); 
        data_dic{mm} = data_dic_cur.DIC3DPPresults; 
        clear data_dic_cur
        %robot data
        FreqFilter_robotdata = 20; 
        if ~isempty(file_robot)
            S = dir(file_robot);
            data_robot{mm} = smp_import(fullfile(S.folder,S.name),1,FreqFilter_robotdata);
        end
        mm = mm+1; 
        elapsedTime = toc;
        fprintf('done in %1.2fsec\n',elapsedTime);
    catch 
        TRIAL_TARGET(mm) = [];  
        data_dic{mm} = [];
        data_robot{mm} = []; 
        fprintf('wrong trial - supressed\n',elapsedTime);
    end
end
Ntrial = length(TRIAL_TARGET); 

fprintf("done\n"); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing & Plot - Force signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
showLogic = 1; 
fprintf('Processing force signals (show = %d)\n',showLogic);

threshold_touch = 0.03; %N
buffer_frame_before = 10; %nbr frame
Nframe_adj = 75; %nbr frame

Nframe = length(data_dic{1}.Points3D);
t = 0:1/actual_FS_vid:(Nframe-1)/actual_FS_vid;
nf_all = zeros(Nframe,Ntrial); 
tf_all = zeros(Nframe,Ntrial); 
% nf_raw = zeros((150+320)*FS_robot/actual_FS_vid,Ntrial); 
% tf_raw = zeros((150+320)*FS_robot/actual_FS_vid,Ntrial); 
nf_raw = zeros((150)*FS_robot/actual_FS_vid,Ntrial); 
tf_raw = zeros((150)*FS_robot/actual_FS_vid,Ntrial); 
idx_touch_all = zeros(1,Ntrial); 

%Init output 
for mm = 1:Ntrial
    data_dic_cur = data_robot{mm};
    %time vector 
    if str2double(cell2mat(regexp(SUBJECT, '\d+', 'match')))<8
        timestart = 6.5;%s
    else
        timestart = 6;%s
    end
    frame2robotMapping = FS_robot/actual_FS_vid; 
    timeStampImg = 1+timestart*FS_robot+((1:Nframe)-1)*frame2robotMapping;
    timeStampRobot = 1+timestart*FS_robot+((1:Nframe*frame2robotMapping)-1);
    nf_cur = data_robot{mm}.nf(timeStampImg);
    tf_cur = data_robot{mm}.tfx(timeStampImg);
    nf_rawcur = data_robot{mm}.nf(timeStampRobot(1)+(1:length(nf_raw)));
    tf_rawcur = data_robot{mm}.tf(timeStampRobot(1)+(1:length(nf_raw)));
    
    %supression of offset 
    nbrimage_baseline = 10; 
    nf_offset = mean(data_robot{mm}.nf(timeStampRobot(1:nbrimage_baseline*frame2robotMapping))); 
    nf_cur = nf_cur - nf_offset;
    nf_rawcur = nf_rawcur - nf_offset;
    
    tf_offset = mean(data_robot{mm}.tfx(timeStampRobot(1:nbrimage_baseline*frame2robotMapping))); 
    tf_cur = tf_cur - tf_offset;
    tf_rawcur = tf_rawcur - tf_offset;
    
    %first instant of loading
    idx_touch = find(nf_cur>threshold_touch,1); 
    
    %output variable 
    nf_all(:,mm) = nf_cur;
    tf_all(:,mm) = tf_cur;
    nf_raw(:,mm) = nf_rawcur;
    tf_raw(:,mm) = tf_rawcur;
    idx_touch_all(mm) = idx_touch;
end

% sync data 
idx_touch_min = min(idx_touch_all);
idx_start = idx_touch_min-buffer_frame_before;
idx_end = idx_start+Nframe_adj;

if showLogic 
% Figure and subplot initialization 
Nrow = 1; 
Ncol = 2;
mergeElement = {};

title_fig = sprintf("%s : Normal force signals %s",SUBJECT,mat2str(TRIAL_TARGET)); 
hf = newfig(title_fig);

set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [ 0    0.1000    0.4000    0.6000] ) ;
 
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 

xvec_frame = 1:Nframe; 
for mm = 1:Ntrial
    nf_cur = nf_all(:,mm);
    idx_touch_cur = idx_touch_all(mm);
    mask_range_cur = idx_touch_cur:idx_touch_cur+Nframe_adj;
    
    ax_idx = 1; 
%     plot(h_ax(ax_idx),data_robot{mm}.time-timestart,data_robot{mm}.nf,'k-');
    plot(h_ax(ax_idx),xvec_frame,nf_cur,'k.-','MarkerSize',10);
    plot(h_ax(ax_idx),xvec_frame(idx_touch_cur),nf_cur(idx_touch_cur),'kx'); 
    title(h_ax(ax_idx),sprintf('trial #%s',num2str(TRIAL_TARGET)));
    h_ax(ax_idx).XLim = [0,Nframe]; 
    h_ax(ax_idx).YLim = [0,6]; 
    
    ax_idx = 2; 
    plot(h_ax(ax_idx),t(mask_range_cur)-t(mask_range_cur(1)),nf_cur(mask_range_cur),'k.','MarkerSize',10);
    if mm == 1; h_ax(ax_idx).YLabel.String = 'Force (N)'; end
    if mm == round(Ntrial/2); h_ax(ax_idx).XLabel.String = 'Time (s)'; end
    title(h_ax(ax_idx),sprintf('Adjusted with %d frames',Nframe_adj));
    h_ax(ax_idx).YLim = [0,6]; 
    
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','on');
grid(h_ax,'on'); 

% Legend
ax_idx = Ncol;
legend_data = cell(1,2);
qw = cell(1,2);
qw{1} = plot(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],'k-');
legend_data{1} = sprintf("%s",'raw signal');
qw{2} = plot(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],'r.');
legend_data{2} = sprintf("%s",'Frame sync');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','southeast');
hleg.Box = 'off';
hleg.NumColumns = 1;
set(hleg,'FontSize',12); 
end


%--------------------------------------------------------------------------
%% Processing - detecting outliers in point cloud 
%--------------------------------------------------------------------------
showLogic = 1; 
fprintf('Processing outliers (show = %d)\n',showLogic);

%Outliers detection on individual maps and 
%Combine individual outliers id. 
mask_outliers_comb = cell(Ntrial,1);
mask_outliers_indiv = cell(Ntrial,1);
id_all = zeros(Ntrial,3); 
for mm = 1:Ntrial
    data_dic_cur = data_dic{mm};
    id = [1,idx_touch_all(mm)+10,length(data_dic_cur.FaceCentroids)];
    
    [mask_comb_current,mask_outliers_indiv_current] = point_cloud_outliers_removal(data_dic_cur,id);
    id_all(mm,:) = id; 
    % Mask outliers combined over multiple instant
    mask_outliers_indiv{mm} = mask_outliers_indiv_current;
    mask_outliers_comb{mm} = mask_comb_current;
    
end

if showLogic
%--------------------------------------------------------------------------
% Figure and subplot initialization 
Nrow = length(id); 
Ncol = Ntrial;
mergeElement = {};

title_fig = sprintf("%s : Map point cloud : outlier detection\n %s",SUBJECT,mat2str(TRIAL_TARGET)); 
hf = newfig(title_fig);

set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.0,0.1,0.4,0.6] ) ;
 
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 

for mm = 1:Ntrial
    data_dic_cur = data_dic{mm};
    id = id_all(mm,:); 
    mask_comb_current = mask_outliers_comb{mm}; 
    for iframeMap = 1:length(id)
        ptFC_cur = data_dic_cur.FaceCentroids{id(iframeMap)};
        mask_indiv_current = mask_outliers_indiv{mm}{iframeMap}; 
        
        ax_idx = mm+(iframeMap-1)*Ncol;
        scatter3(h_ax(ax_idx),ptFC_cur(:,1),...
            ptFC_cur(:,2),...
            ptFC_cur(:,3),...
            1,'Marker','.',...
            'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
        scatter3(h_ax(ax_idx),ptFC_cur(mask_indiv_current,1),...
            ptFC_cur(mask_indiv_current,2),...
            ptFC_cur(mask_indiv_current,3),...
            20,'Marker','.',...
            'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter3(h_ax(ax_idx),ptFC_cur(mask_comb_current&~mask_indiv_current,1),...
            ptFC_cur(mask_comb_current&~mask_indiv_current,2),...
            ptFC_cur(mask_comb_current&~mask_indiv_current,3),...
            20,'Marker','.',...
            'MarkerEdgeColor','b','MarkerFaceColor','b');
        h_ax(ax_idx).DataAspectRatio = [1 1 1];
        
        % Axis properties
        if mm == 1
            text(h_ax(ax_idx),min(ptFC_cur(:,1),[],1,'omitnan'),mean(ptFC_cur(:,2),1,'omitnan'),0,sprintf("frame #%d",id(iframeMap)),...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',15,...
                'FontWeight','bold');
        end
        if iframeMap==1; title(h_ax(ax_idx),sprintf('trial #%d',TRIAL_TARGET(mm))); end
        xlim(h_ax(ax_idx),[-10,10]);
        ylim(h_ax(ax_idx),[-5,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','on');
        set(findall(h_ax(ax_idx), 'type', 'text'), 'visible', 'on')
    end
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing - data collection, masking and filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Processing - data collection, masking and filtering\n');
optStruct = struct(); 
optStruct.corr_threshold = 1; 
optStruct.freqTempFilt = 10; 
optStruct.smoothPar.n=4; % SPACE FILTERING PARAM 1
optStruct.smoothPar.sigma=100; % SPACE FILTERING PARAM 2
Nface = length(FACEMEASURE);

T = cell(Ntrial,1); 

for mm = 1:Ntrial
    tic; 
    fprintf('--> trial %d...',TRIAL_TARGET(mm));
    data_dic_cur = data_dic{mm};
    FC_all = struct(); 
    
    param = []; 
    % Retreive necessary data
    Nframe = length(data_dic_cur.FaceCentroids);
    Fcentroidfirst=data_dic_cur.FaceCentroids{1};
    Fcentroidlast=data_dic_cur.FaceCentroids{end};
    Npoints = size(Fcentroidfirst,1);
    
    %get smooth space parameter for the first pc map 
    [~,idx_nearest_neighbours_set,weight_nearest_neighbours_set]=smooth_space(false(Npoints,1),Fcentroidfirst,optStruct);
    
    % mask badly correlated faces  
    corrLogic = data_dic_cur.FaceCorrComb{end} > optStruct.corr_threshold; 
    % mask stitch faces
    gapLogicStitch = data_dic_cur.FacePairInds==3;%|r.FacePairInds==1;% ATTENTION A CHANGER !!
    % mask closest points to stitch faces 
    gapLogicnearest = false(Npoints,1); 
    gapLogicnearest(idx_nearest_neighbours_set(gapLogicStitch,2)) = true; %first closest neighbor 
    gapLogicnearest(idx_nearest_neighbours_set(gapLogicStitch,3)) = true; %second closest neighbor 
    % mask outliers 
    outliersLogic = mask_outliers_comb{mm};
    % mask combined of point to suppress
    mask_pc2suppress = corrLogic|gapLogicStitch|gapLogicnearest|outliersLogic;
    Fcentroidfirst(mask_pc2suppress,:)=NaN; 
    Fcentroidlast(mask_pc2suppress,:)=NaN;  
    
    % time filtering 
    data_dic_cur = filter_measures(data_dic_cur,optStruct.freqTempFilt,actual_FS_vid,FACEMEASURE(1:length(FACEMEASURE)-1)); %last is curvature
        
    for is = 1:Nface-1 % without curvature measurements here.
        faceMeasure_is = FACEMEASURE{is}; 
        
        %Load current measurements from DIC analysis
        [FCmat_cum,FCmat_rate] = data_dic2matrix(data_dic_cur,faceMeasure_is); 
        if strcmp(faceMeasure_is,'Area')
            A = FCmat_cum(:,1);
            [FCmat_cum,FCmat_rate] = data_dic2matrix(data_dic_cur,'J'); 
            FCmat_cum = A.*(1+FCmat_cum);
            FCmat_rate = A.*(1+FCmat_rate);
        elseif strcmp(faceMeasure_is,'Epc1')||strcmp(faceMeasure_is,'Epc2')
            FCmat_cum = FCmat_cum*1e2; %percent
            FCmat_rate = FCmat_rate*1e2*actual_FS_vid; %percent per sec
            
            %Collecting principal directions
            FC_all.(faceMeasure_is).vec = data_dic_cur.Deform.(sprintf('%svecCur',faceMeasure_is)){Nframe};
            FC_all.(faceMeasure_is).vec(mask_pc2suppress,:) = NaN;
        end
        
        FCmat_cum(mask_pc2suppress,:) = NaN; 
        FCmat_rate(mask_pc2suppress,:) = NaN; 
        
        %gather data of each faces
        FC_all.(faceMeasure_is).(deftypeNames{1}) = FCmat_cum;
        FC_all.(faceMeasure_is).(deftypeNames{2}) = FCmat_rate;
        
        % loop on every deformation type and on everyframe
        for iframe = 1:Nframe
            Fcentroidnow=data_dic_cur.FaceCentroids{iframe}; %Npoints x 3
            Fcentroidnow(mask_pc2suppress,:) = NaN; 
            FC_all.xyz(:,:,iframe) = Fcentroidnow;
            for idef = 1:length(deftypeNames)
                % Retreive current data
                FCnow = FC_all.(faceMeasure_is).(deftypeNames{idef})(:,iframe);
                % Face spatial filtering
                FCnow=smooth_space(FCnow,Fcentroidnow,optStruct,...
                    'idx_nearest_neighbours_set',idx_nearest_neighbours_set,...
                    'weight_nearest_neighbours_set',weight_nearest_neighbours_set);
                % output variable
                FC_all.(faceMeasure_is).(deftypeNames{idef})(:,iframe)=FCnow;
            end
        end
    end
    param.Nsuppressed = sum(mask_pc2suppress);
    param.idx_nearest_neighbours_set = idx_nearest_neighbours_set;
    param.weight_nearest_neighbours_set = weight_nearest_neighbours_set; 

    %Gather data into structure
    T{mm}.trial_nbr = TRIAL_TARGET(mm); 
    T{mm}.FRcond = TRIAL_TARGET(mm); 
    T{mm}.Nframe = Nframe; 
    T{mm}.nfcond = trial_info.Nfcond(TRIAL_TARGET(mm) == trial_info.TrialNbr);
    T{mm}.spddxlcond = trial_info.Spddxlcond(TRIAL_TARGET(mm) == trial_info.TrialNbr);
    T{mm}.FC = FC_all;
    T{mm}.nf = nf_all(:,mm); 
    T{mm}.nf_touch_idx = idx_touch_all(mm); 
    T{mm}.tf = tf_all(:,mm); 
    T{mm}.param.specific = param; 
    T{mm}.param.common = optStruct; 
    T{mm}.mask_supress = mask_pc2suppress; 
    T{mm}.Npoints = size(Fcentroidlast,1); 

    elapsedTime = toc;
    fprintf('%1.2fsec\n',elapsedTime);
end

%Clear unnecessary output
clear FC_all r 

fprintf('done\n');

%% Processing - Segmentation and contact region
% REM : We add curvature measurement here. 
fprintf('Processing - segmentation and contact region\n');

optStruct.Nregion = 5; 
optStruct.guess_threshold = -1; %mm 
optStruct.min_inside = 3; %#point inside 
optStruct.margin_around_plane = 0.2; %mm 
optStruct.border_width_annulus = 1; %mm
optStruct.connex_factor = 0.1; % connex factor for the boundary location of contact border
optStruct.numNeighbours = 20; 
optStruct.threshold_percent_last = 75; 

for mm = 1:Ntrial
    tic; 
    fprintf('--> trial %d...',TRIAL_TARGET(mm));
    Nframe = T{mm}.Nframe; 
    Fcentroidlast = T{mm}.FC.xyz(:,:,end);
    Npoints = T{mm}.Npoints;
    curvature_all = zeros(Npoints,Nframe);
    
    %computation of z-threshold from last map point cloud
    temp = fit_xyplane_below_threshold_pc(Fcentroidlast,optStruct.guess_threshold,'margin',optStruct.margin_around_plane);
    param.zthreshold_lastmap = mean(Fcentroidlast(temp,3));
    
    %computation of absolute threshold for last map curvature
    temp = Fcentroidlast;
    temp(Fcentroidlast>param.zthreshold_lastmap) = NaN; 
    [~,param.threshold_absolute_set] = find_curv_threshold(Fcentroidlast(Fcentroidlast(:,3)<param.zthreshold_lastmap,:),'numNeighbours',optStruct.numNeighbours,'threshold_percent',optStruct.threshold_percent_last);

    % data allocation
    mask_false = false(Npoints,Nframe);
    mask_in_all = mask_false;
    mask_in_below_plane_all = mask_false;
    mask_annulus_all = mask_false;
    mask_border_all = mask_false;
    CA_all = cell(1,Nframe);
    for iframe = 1:Nframe 
        Fcentroidnow = T{mm}.FC.xyz(:,:,iframe);
        
        % Map current curvature measurement : based on geometric
        % consideration of the object, its point cloud local variance 
        % Rem : the mask has already been applied to the xyz point cloud
        % before.
        [mask_curv_now,threshold_now,curvature_now] = find_curv_threshold(Fcentroidnow,...
            'numNeighbours',optStruct.numNeighbours,...
            'threshold_absolute',param.threshold_absolute_set,...
            'idx_target_curv',1);
    
        % Contact region
        [mask_in_now,mask_in_below_plane_now] = find_contact_area(Fcentroidnow,Fcentroidlast,mask_curv_now,...
            'zthreshold_lastmap',param.zthreshold_lastmap,...
            'margin_around_plane',optStruct.margin_around_plane,...
            'threshold_absolute_set',param.threshold_absolute_set);

        % Border of contact 
        CA_now = [NaN,NaN,NaN];
        mask_border_now = false(Npoints,1);
        if sum(mask_in_now) > optStruct.min_inside
            temp = Fcentroidnow; temp(~mask_in_now,:) = NaN; 
            [CA_now,mask_border_now] = findboundary(temp,'connex_factor',optStruct.connex_factor);
        end
        % Annulus around contact border
        mask_in_annulus_now = find_annulus_border(Fcentroidnow,CA_now,'borderWidth',optStruct.border_width_annulus);

        %Last frame computation
        if iframe == Nframe
            % Segmentation of the fingertip considering the
            % same part for the entire trial.
            % Meaning that each segmentation represent an
            % actual part of the fingertip that we follow.
            [mask_region,mask_RIN1,mask_RIN2,mask_RIN3,mask_out] ...
                = find_fingerside(Fcentroidfirst,Fcentroidlast,CA_now,optStruct.Nregion);
            CA_last = CA_now;
        end
        % Curvature spatial filtering
        curvature_now=smooth_space(curvature_now,Fcentroidnow,optStruct,...
            'idx_nearest_neighbours_set',T{mm}.param.specific.idx_nearest_neighbours_set,...
            'weight_nearest_neighbours_set',T{mm}.param.specific.weight_nearest_neighbours_set);
        
        mask_in_all(:,iframe) = mask_in_now; 
        mask_in_below_plane_all(:,iframe) = mask_in_below_plane_now;
        mask_annulus_all(:,iframe) = mask_in_annulus_now; 
        mask_border_all(:,iframe) = mask_border_now;
        curvature_all(:,iframe) = curvature_now;
        CA_all{iframe} = CA_now;
    end
    
    T{mm}.param.specific.zthreshold_lastmap = param.zthreshold_lastmap; 
    T{mm}.param.specific.threshold_absolute_set = param.threshold_absolute_set; 
    T{mm}.FC.curv.cum = curvature_all;
    T{mm}.FC.curv.rate = curvature_all;
        
    % Contact region circumpherencial
    Ncirc = size(mask_region,2); %Nbr segmentatin radial
    mask_R_all = cell(1,Ncirc);
    for ireg = 1:Ncirc
        mask_R_all{ireg} = repelem(mask_region(:,ireg),1,Nframe);
    end
    
    % Contact segmentation : first in/early in/last in/out
    % Number of segmentation
    mask_all = true(size(mask_RIN1));
    mask_ROUT = repelem(mask_all&~mask_RIN3&~mask_RIN2&~mask_RIN1,1,Nframe);
    mask_RIN1 = repelem(mask_RIN1,1,Nframe);
    mask_RIN2 = repelem(mask_RIN2,1,Nframe);
    mask_RIN3 = repelem(mask_RIN3,1,Nframe);
    mask_out = repelem(mask_out,1,Nframe);
    
    T{mm}.name_circ = {'west','north-west','north','north-east','east'};
    T{mm}.name_rad = {'first in','early_in','last in','outside','extreme_out'};
    T{mm}.mask_circ = mask_R_all;%{mask_R1_all,mask_R2_all,mask_R3_all,mask_R4_all,mask_R5_all};
    T{mm}.mask_rad = {mask_RIN1,mask_RIN2,mask_RIN3,mask_ROUT,mask_out}; 
    T{mm}.mask_in = mask_in_all;
    T{mm}.mask_in_below_plane = mask_in_below_plane_all;
    T{mm}.mask_in_annulus = mask_annulus_all;
    T{mm}.mask_in_border = mask_border_all;
    T{mm}.ca = CA_all; 
    
    elapsedTime = toc;
    fprintf('%1.2fsec\n',elapsedTime);
end
clear mask_R_all mask_all mask_ROUT mask_RIN1 mask_RIN2 mask_RIN3 mask_R0_all
Nrad = length(T{mm}.mask_rad); 

fprintf('done\n');

%% Angle of attack measurement 
clim_set = [0,0.01];%0.05; 
list_force_target = [0.05,0.5,1:0.5:5]; 
gridStep = 1;%0.75;
numNeighbours = 8; 
scatter_width = 20;

angle_saved_plot = NaN(Ntrial,length(list_force_target));
showLogic = true; 

if showLogic 
h_fig = newfig(sprintf('Attack angle measurement %s',SUBJECT)); 
set(h_fig, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.0005    0.7444    0.9979    0.1769]);
Nrow = length(list_force_target); Ncol = Ntrial; 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
end

for mm = 1:Ntrial
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>NFCOND_SET,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    
    id = list_frame_num; 
    
    PCcurrent = T{mm}.FC.xyz; 
    axis_target = 1; 
    ii = 0; 
    for iframe = id
        ii = ii+1;
        
        ptFC_cur = PCcurrent(:,:,iframe);
        
        logicNaN = isnan(ptFC_cur(:,1));
        points = ptFC_cur(~logicNaN,:);
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

        %1)take the boundary
        %2)fine interpolation 
        %3)fixed inter-distance between points 
        k = boundary(pointsds(:,other_axis(1)),pointsds(:,other_axis(2)),0.1);
        pointsds = pointsds(k,:);
        [pointsds,dudt,fofthandle] = interparc(100,pointsds(:,1),pointsds(:,2),pointsds(:,3),'linear');
        temp = pointCloud(pointsds);
        temp = pcdownsample(temp,'gridAverage',gridStep);
        pointsds = temp.Location;

        %find the curvature
        [~,~,curvature_now] = find_curv_threshold(pointsds,...
            'numNeighbours',numNeighbours,...
            'threshold_absolute',0,...
            'idx_target_curv',2);
        curv_lim_set = 0.0005;
        zlim_set = 2; 
        ylim_set = 5; 
        xvec = -10:25; 
        mask_curv = curvature_now < curv_lim_set; 
        y = pointsds(:,2); 
        z = pointsds(:,3); 
        maskgeom1 = z>zlim_set & y<ylim_set; 
        maskgeom2 = z<zlim_set & y<ylim_set; 
        mask1 = mask_curv & maskgeom1; 
        mask2 = mask_curv & maskgeom2; 

        %mask 1
        x = pointsds(mask1,2); y = pointsds(mask1,3);
        c1 = polyfit(x,y,1); y_est1 = polyval(c1,xvec);
        %mask 2
        x = pointsds(mask2,2); y = pointsds(mask2,3);
        c2 = polyfit(x,y,1); y_est2 = polyval(c2,xvec);
        
        angle = (atan(c1(1))-atan(c2(1)))*180/pi; 
%         angle = (atan(c1(1)))*180/pi; 
        angle_saved_plot(mm,ii) = angle;
        if ii == 1
            % Gather information 
            T{mm}.angle = angle; 
        end

        if showLogic
        % graph 
        ax_idx = mm+(ii-1)*Ncol; 
        %map 
        scatter3(h_ax(ax_idx),points(:,1),points(:,2),points(:,3),1,'filled');
        scatter3(h_ax(ax_idx),pointsds(:,1),pointsds(:,2),pointsds(:,3),scatter_width,curvature_now,'filled');
        %line fit 
        plot3(h_ax(ax_idx),zeros(length(xvec),1),xvec,y_est1,'k');
        h_line = plot3(h_ax(ax_idx),zeros(length(xvec),1),xvec,zeros(length(xvec),1),'k');
        
        text(h_ax(ax_idx),15,0,15,sprintf('trial #%d\nangle = %1.2f°\nframe %d',TRIAL_TARGET(mm),angle,iframe)); 
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
        set(h_ax(ax_idx),'clim',clim_set)
        xlim(h_ax(ax_idx),[-15,15]);
        ylim(h_ax(ax_idx),[-10,25]);
        zlim(h_ax(ax_idx),[-3,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View',view_plot);
        colormap(h_ax(ax_idx),parula);
%         colorbar(h_ax(ax_idx)); 
        end
    end
end


if showLogic
h_fig = newfig(sprintf('Attack angle measurement evolution%s',SUBJECT)); 
set(h_fig, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.0005    0.7444    0.9979    0.1769]);
Nrow = 1; Ncol = 1; 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 

for mm = 1:Ntrial
    ax_idx = 1;
    plot(h_ax(ax_idx),list_force_target,angle_saved_plot(mm,:),'.-');
end
set(h_ax(ax_idx),'YLim',[0,25]);
end
fprintf('done\n');

%% SPEED condition 

color_spd = {[1 1 1].*0.5,[1 1 1].*0.0};

logicspd = zeros(1,Ntrial); 
P = zeros(2,Ntrial);

fig_name = sprintf("%s - Speed condition",SUBJECT); 
h_fig = newfig(fig_name); set(h_fig,'Unit','Normalized','Position',[0.6901    0.6111    0.2458    0.2796]);
h_ax = subplot_ax(1,2); hold(h_ax,'on'); grid(h_ax,'on'); set(h_ax,'FontSize',12); 
for mm = 1:Ntrial
    spdcond = T{mm}.spddxlcond;
    logicspd(mm) = spdcond==0.08; 
    
    frame_vec = 1:T{mm}.Nframe;
    idx_s = T{mm}.nf_touch_idx;
    nf_cur = T{mm}.nf(idx_s:end);
    idx1 = find(nf_cur >= 0.1*T{mm}.nfcond,1,'first');
    idx2 = find(nf_cur >= 0.9*T{mm}.nfcond,1,'first');
    nf_cur = nf_cur(1:idx2+50);
    
    frame_span = idx1:idx2;
    P(:,mm) = polyfit(frame_span,nf_cur(frame_span),1);
    
    T{mm}.nfslope = P(1,mm)*50; 
    
    yfit = P(1,mm)*frame_span+P(2,mm);
    ax_idx = 1; 
    hp = plot(h_ax(ax_idx),(1:length(nf_cur))./actual_FS_vid,nf_cur,'Color',color_spd{logicspd(mm)+1});
    plot(h_ax(ax_idx),[idx1,idx2]./actual_FS_vid,nf_cur([idx1,idx2]),'k.','MarkerSize',10);
    plot(h_ax(ax_idx),frame_span./actual_FS_vid,yfit,'-','Color',hp.Color,'LineWidth',2);
    
end
ax_idx = 1; 
qw = cell(1,2); 
for ispd = 1:length(color_spd)
    qw{ispd} = plot(h_ax(ax_idx),NaN,NaN,'Color',color_spd{ispd},'LineWidth',2);
end
h_leg = legend(h_ax(ax_idx),[qw{:}],["low","high"],'Position',[0.2915,0.1947,0.1674,0.1407],'box','off'); 
h_ax(ax_idx).XLim = [0,1.50]; 
xlabel(h_ax(ax_idx),'Time (s)'); 
ylabel(h_ax(ax_idx),'NF (N)'); 

ax_idx = 2; 
scatter(h_ax(ax_idx),logicspd,P(1,:)*50,20,hp.Color,"filled",'XJitter','randn','XJitterWidth',0.3); 
h_ax(ax_idx).XTick = [0,1]; 
h_ax(ax_idx).XTickLabel = {'low','high'}; 
h_ax(ax_idx).XLim = [-0.50,1.50]; 
h_ax(ax_idx).YLim = [0,15]; 
title(h_ax(ax_idx),["Differenciation between","speed condition"]); 
xlabel(h_ax(ax_idx),'Speed condition'); 
ylabel(h_ax(ax_idx),'Fitted slope (N/s)'); 

exportgraphics(h_fig,sprintf('fig/%s.%s',fig_name,'png')); 
fprintf('done\n'); 

%% Save current data
fprintf('Saving data...[');
for mm = 1:Ntrial
    trial_cur = T{mm}.trial_nbr; 
    folder_loc = fullfile(baseResultPath,'analysis',SUBJECT); 
    if ~exist(folder_loc,'dir')
        mkdir(folder_loc);
    end
    filename = sprintf('%s_%03d_proccessed_data_loading_%s.mat',SUBJECT,trial_cur,FRcond{MATERIAL_NBR}); 
    data = T{mm};
    save(fullfile(folder_loc,filename),'data');
    fprintf('%d ',trial_cur);
end
fprintf('] '); 
fprintf('done\n');

%% Load processed data
fprintf('Load data...[');
T = cell(Ntrial,1); 
itrial = 0; 
for mm = 1:Ntrial
    trial_cur = TRIAL_TARGET(mm);
    filename = sprintf('%s_%03d_proccessed_data_loading_%s.mat',SUBJECT,trial_cur,FRcond{MATERIAL_NBR});
    try
        data = load(fullfile(baseResultPath,'analysis',SUBJECT,filename)); 
        itrial = itrial+1; 
        T{itrial} = data.data;
        fprintf('%d ',trial_cur);
    catch
    end
end
Ntrial = itrial; 
fprintf('] '); 
fprintf('done\n');

%% PLOT ALL FRAME OF A SINGLE MEASUREMENT MAP 
face_target = 4; 
trial_target = 1;
deftype = 'rate';

for mm = trial_target
    for iface = face_target 
        faceMeasure_ii = FACEMEASURE{iface}; 
        id = 1:57;%size(T{mm}.FC.xyz,3);
        xyz = T{mm}.FC.xyz;
        mask_in = T{mm}.mask_in; 
        val = T{mm}.FC.(faceMeasure_ii).(deftype);
        title_fig = sprintf('%s#%03d %s %s all filt',SUBJECT,T{mm}.trial_nbr,deftype,faceMeasure_ii); 
        h_fig = newfig(title_fig); hold on;
        set(h_fig,'Position',[0    0.7938   50.8000   25.7440]); 
        Nrow = 3;%8; 
        Ncol = 19; 
        merge_element = {};
        ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(ax,'on'); 
        for iframe = id
            ax_idx = iframe; 
            
            scatter3(ax(ax_idx),...
                xyz(:,1,id(iframe)),...
                xyz(:,2,id(iframe)),...
                xyz(:,3,id(iframe)),...
                2,...
                val(:,id(iframe)),...
                'filled');
            plot3(ax(ax_idx),...
                T{mm}.ca{iframe}(:,1),...
                T{mm}.ca{iframe}(:,2),...
                T{mm}.ca{iframe}(:,3),...
                'k-','LineWidth',2); 
            
            % Border of contact 
            mask_in_now = T{mm}.mask_in_below_plane(:,id(iframe)); 
            CA_now = T{mm}.ca{iframe};
            plot3(ax(ax_idx),...
                CA_now(:,1),...
                CA_now(:,2),...
                CA_now(:,3),...
                '-','LineWidth',1,'Color',[1 1 1].*0.5); 
                
            text(ax(ax_idx),min(xyz(:,1,id(iframe)),[],1,'omitnan'),mean(xyz(:,2,id(iframe)),1,'omitnan'),0,sprintf('%03d',id(iframe)),...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',12);
            xlim(ax(ax_idx),[-10,10]);
            ylim(ax(ax_idx),[-5,15]);
            colormap(ax(ax_idx),'jet');
            set(ax(ax_idx),'clim',[0,2e-3]); %,[-0.10,0.10]); 
            set(ax(ax_idx),'FontSize',12,'clipping','off');
            set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90],'visible','off'); 
        end
        for ii = length(id)+1:length(ax)
            ax_idx = ii; 
            set(ax(ax_idx),'visible','off'); 
        end
    end
end

%% Alternative/complementary measurement of contact area using displacement rate 
nepochs = 20; 
start_looking = 25; 
for mm = 1
% title_fig = sprintf('Search for the first contact moment trial %d',T{mm}.trial_nbr); 
% newfig(title_fig); hold on;
% ax = subplot_ax(2,1); hold(ax,'on'); 
    % Normalized threshold of the ratio between max and min of rate
    % displacement : Define how far from the rest of the point, the minimum
    % displacement is from the current displacement of the bunch. 
    norm_disp_rate_threhold_pc = 0.9; 
    
    xyz = T{mm}.FC.xyz;
    val = abs(T{mm}.FC.DispZ.rate);
    
    % Point sub-guess-threshold of -1
    near_contact_logic = squeeze(xyz(:,3,:)) < -0.50; 
    start_looking = find(sum(near_contact_logic,1)>1,1); 
    
    % Masking other points
    val(~near_contact_logic) = NaN;
    % 
    meanval = prctile(val,50,1);
    val_current_norm = val./meanval; 
    disp_rate_low = prctile(val_current_norm,0.1,1);
    logic_contact = disp_rate_low < norm_disp_rate_threhold_pc;
    
    
%     %ax(1)
%     ax_idx = 1; 
%     hp=plot(ax(ax_idx), [false,logic_contact(2:end)]); 
%     xline(ax(ax_idx),start_looking,'-','Color',get(hp,'Color'));
%     
%     %ax(2)
%     ax_idx = 2; 
%     plot(ax(ax_idx), disp_rate_low,'Color',get(hp,'Color'))
%    
%     yline(ax(ax_idx), norm_disp_rate_threhold_pc);  
%     linkaxes(ax,'x'); 

    
    val = T{mm}.FC.DispZ.rate;
    
    clim_set = [-0.129,0.129]; 
    id = start_looking:start_looking+nepochs-1;
    
    title_fig = sprintf('Search for the first contact map trial %d',T{mm}.trial_nbr); 
    h_fig = newfig(title_fig); hold on;
    set(h_fig,'Position',[5.3181   11.9856   25.6381   12.5942]); 
    Nrow = 2; Ncol = length(id); 
    merge_element = {1:Ncol};
    ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(ax,'on'); 
    for iframe = 1:length(id)
        near_contact_logic_current = near_contact_logic(:,id(iframe)); 
        val_current = val(:,id(iframe)); 
        logic_contact_current = logic_contact(id(iframe)); 
        logic_disp_rate_low_points = abs(val_current) <= repelem(0.04,size(val_current,1),1); 
        ax_idx = iframe+0*Ncol; 
        if iframe == 1
            timeframe = start_looking:start_looking+nepochs+5;
            yline(ax(ax_idx), norm_disp_rate_threhold_pc); 
            plot(ax(ax_idx),timeframe,disp_rate_low(timeframe),'k.-','MarkerSize',10,'LineWidth',2);
            xlim(ax(ax_idx),[timeframe(1),timeframe(end)]); 
            ax(ax_idx).XTickLabel = [];
            set(get(ax(ax_idx), 'XAxis'), 'Visible', 'off');
            ylabel(ax(ax_idx),'Disparity dispZ rate'); 
            ax(ax_idx).YGrid = 'on';
            ax(ax_idx).GridLineStyle = '-';
            text(ax(ax_idx),find(logic_contact,1)+0.2,disp_rate_low(find(logic_contact,1)),sprintf('#%d',find(logic_contact,1)),...
                'Rotation',0.,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',12);
        end
        plot(ax(ax_idx),[id(iframe),id(iframe)],[0,disp_rate_low(id(iframe))],'k-','LineWidth',1);
        if ~logic_contact_current
            plot(ax(ax_idx),id(iframe),disp_rate_low(id(iframe)),'ks','MarkerSize',10,'MarkerFaceColor',[1 1 1].*0.8);
        else
            plot(ax(ax_idx),id(iframe),disp_rate_low(id(iframe)),'ks','MarkerSize',10,'MarkerFaceColor','k');
        end
   
        ax_idx = iframe+1*Ncol; 
%         title(ax(ax_idx),sprintf('#%d',id(iframe))); 
        scatter(ax(ax_idx),xyz(:,1,id(iframe)),...
            xyz(:,2,id(iframe))+30,...
            1,'filled'); 
        scatter(ax(ax_idx),xyz(logic_disp_rate_low_points,1,id(iframe)),...
            xyz(logic_disp_rate_low_points,2,id(iframe))+30,...
            2,'r','filled'); 
        
        scatter(ax(ax_idx),xyz(:,1,id(iframe)),...
            xyz(:,2,id(iframe)),...
            2,val_current,'filled'); 
        if iframe == 1
            text(ax(ax_idx),min(xyz(:,1,id(iframe)),[],1,'omitnan'),mean(xyz(:,2,id(iframe)),1,'omitnan'),0,'d/dt Disp Z',...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',12);
            text(ax(ax_idx),min(xyz(:,1,id(iframe)),[],1,'omitnan'),mean(xyz(:,2,id(iframe)),1,'omitnan')+30,0,'Below',...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',12);
        end
%         x_axtrace = id(iframe); y_axtrace = 0;
%         x_axmap = mean(xyz(:,1,id(iframe)),1,'omitnan'); y_axmap = min(xyz(:,2,id(iframe)),[],1,'omitnan'); 
%         h = line([x_axmap x_axtrace],[y_axmap y_axtrace]);
%         set(h,'LineWidth',2,'Color','k')
        colormap jet 
        set(ax(ax_idx),'clim',clim_set); 
        xlim(ax(ax_idx),[-10,10]);
        ylim(ax(ax_idx),[-5,15+30]);
%         zlim(ax(ax_idx),[-2,15]);%
        set(ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90],'visible','off'); 
        set(findall(ax(ax_idx), 'type', 'text'), 'visible', 'on')
    end
    set(ax,'FontSize',12,'clipping','off');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAP - point cloud 
%ax last col° : 3D fingertip segmentation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color_contraction = [1 1 1].*0.5;
color_elongation = [1 1 1].*0.2;
Nrad = length(T{1}.mask_rad); 
Ncirc = length(T{1}.mask_circ); 
%--------------------------------------------------------------------------
% Segmentation RADIAL 
%--------------------------------------------------------------------------

% settings : colormap, legend name,...
C=[0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
color_set_array = C(1:Nrad,:);
marker_set = 'o'; 
lgd_name = T{1}.name_rad;
color_set_cell = mat2cell(color_set_array,repelem(1,size(color_set_array,1)),[3]);

% Figure and subplot initialization 
Nrow = Ncirc; 
Ncol = Ntrial+1;
mergeElement = {};

title_fig = sprintf("%s : Map point cloud segmentation radial\n%s",SUBJECT,mat2str(TRIAL_TARGET)); 
hf = newfig(title_fig);

set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.0,0.1,0.4,0.6] ) ;
 
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 

for mm = 1:Ntrial
Nframe = T{mm}.Nframe; 
id = Nframe; 


for iframeMap = 1:length(id)
    ptFC_cur = T{mm}.FC.xyz(:,:,id(iframeMap)); %r.FaceCentroids{id(iframeMap)};
    FC_epc1 = T{mm}.FC.Epc1.cum(:,id(iframeMap)); D1 = T{mm}.FC.Epc1.vec;
    FC_epc2 = T{mm}.FC.Epc2.cum(:,id(iframeMap)); D2 = T{mm}.FC.Epc2.vec;
    for xx = 1:Ncirc
        ax_idx = mm+(xx-1)*Ncol;
        mask_circ_current = T{mm}.mask_circ{xx};
        scatter(h_ax(ax_idx),ptFC_cur(:,1),...
            ptFC_cur(:,2),...
            1,'Marker','.',...
            'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
        
        
        h_quiv1 = gobjects(1,Nrad);
        h_quiv2 = gobjects(1,Nrad);
        for kk =1:Nrad
            mask_rad_current = T{mm}.mask_rad{kk};
            current_mask = mask_rad_current(:,1)&mask_circ_current(:,1);
            xyz = ptFC_cur(current_mask,:);
            FC_epc1_region = FC_epc1(current_mask,:);
            FC_epc2_region = FC_epc2(current_mask,:);
            D1_region = D1(current_mask,:);
            D2_region = D2(current_mask,:);
            
            % fit ellipse from last contact contour
            [~,centerEl] = myfitellipse(T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2));
            zplane_mean = mean(T{mm}.ca{end}(:,3)); 
            centerEl = [centerEl,zplane_mean];
            pnow = mean(xyz,1,'omitnan'); %centroid

            vecunit_now = pnow-centerEl;
            product_vec1 = dot(repmat(vecunit_now(1:2),length(D1_region),1),D1_region(:,1:2),2);
            product_vec2 = dot(repmat(flip(vecunit_now(1:2)),length(D2_region),1),D2_region(:,1:2),2);
            orientation1 = (product_vec1>0)*2-1; %positive = +1 / negative = -1
            orientation2 = (product_vec2>0)*2-1; %positive = +1 / negative = -1
            D1_region(:,1:2) = D1_region(:,1:2).*orientation1;
            D2_region(:,1:2) = D2_region(:,1:2).*orientation2;
            Dsnow1 = median(D1_region,1,'omitnan').*abs(median(FC_epc1_region,1,'omitnan'))/2; %scaled vector 
            Dsnow2 = median(D2_region,1,'omitnan').*abs(median(FC_epc2_region,1,'omitnan'))/2; %scaled vector                
            
            %plot X-Y points of current region 
            scatter(h_ax(ax_idx),xyz(:,1),xyz(:,2),5,'Marker',marker_set,...
                'MarkerEdgeColor',color_set_cell{kk},'MarkerFaceColor',color_set_cell{kk});
            %plot principales directions projected onto 2D plane 
            h_quiv1(kk) = quiver(h_ax(ax_idx),pnow(1)-Dsnow1(:,1)/2,pnow(2)-Dsnow1(:,2)/2,Dsnow1(:,1),Dsnow1(:,2),0,...
                'Color',color_contraction,'ShowArrowHead','off','AutoScale','off','LineWidth',2);
            h_quiv2(kk) = quiver(h_ax(ax_idx),pnow(1)-Dsnow2(:,1)/2,pnow(2)-Dsnow2(:,2)/2,Dsnow2(:,1),Dsnow2(:,2),0,...
                'Color',color_elongation,'ShowArrowHead','off','AutoScale','off','LineWidth',2);
            
        end
        % plot contour 
        plot(h_ax(ax_idx),T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2),'k-','LineWidth',1,'Color',[1 1 1].*0.3);
        
        % title(h_ax(ax_idx),region_name_set{xx});
        % xlabel(h_ax(ax_idx),'x'); ylabel(h_ax(ax_idx),'y');
        if mm == 1 
            text(h_ax(ax_idx),min(ptFC_cur(:,1),[],1,'omitnan'),mean(ptFC_cur(:,2),1,'omitnan'),0,T{mm}.name_circ{xx},...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',15,...
                'FontWeight','bold');
        end
        if xx == 1
            text(h_ax(ax_idx),mean(ptFC_cur(:,1),1,'omitnan'),min(ptFC_cur(:,2),[],1,'omitnan'),0,sprintf('trial #%d',T{mm}.trial_nbr),...
            'Rotation',0,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15,...
            'FontWeight','bold');
        end
        xlim(h_ax(ax_idx),[-10,10]);
        ylim(h_ax(ax_idx),[-5,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');    
%         set(findall(h_ax(ax_idx), 'type', 'text'), 'visible', 'on')

        % set quiver plot to front of the axis.
        uistack(h_quiv1,'top');
        uistack(h_quiv2,'top');
    end
end
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off','visible','off','Color','None');

% Legend - 1°
ax_idx = 2*Ncol; 
legend_data = cell(1,2);
qw = cell(1,2);
qw{1} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_contraction);
legend_data{1} = sprintf("%s",'Min ~ Contraction');
qw{2} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_elongation);
legend_data{2} = sprintf("%s",'Max ~ Elongation');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','northeast');
hleg.Box = 'off';
hleg.NumColumns = 1;
hleg.Title.String = 'Princ. dir.'; 
hleg.Orientation = 'horizontal';
set(hleg,'FontSize',12); 

% Legend - 2°
ax_idx = Ncol;%2*NplotTimeSize; 
legend_data = cell(1,Nrad+1);
qw = cell(1,Nrad+1);
for kk = 1:Nrad
    qw{kk} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_set_cell{kk});
    legend_data{kk} = sprintf("%s",lgd_name{kk});
end
qw{kk+1} = plot(h_ax(ax_idx),NaN,NaN,'-','Color',[1 1 1].*0.4);
legend_data{kk+1} = sprintf("%s",'Last contour');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','best','Interpreter','none');
hleg.Box = 'off';
hleg.NumColumns = 1;
hleg.Title.String = 'Region'; 
set(hleg,'FontSize',12); 


%--------------------------------------------------------------------------
% Segmentation CIRCUMPHERENCIAL 
%--------------------------------------------------------------------------

% settings : colormap, legend name,...
C = colororder(autumn(Ncirc)); %parula
color_set_array = C(1:Ncirc,:);
marker_set = 'o'; 
lgd_name = T{1}.name_circ;
color_set_cell = mat2cell(color_set_array,repelem(1,size(color_set_array,1)),[3]);

% Figure and subplot initialization 
title_fig = sprintf("%s : Map point cloud segmentation circumpherencial\n%s",SUBJECT,mat2str(TRIAL_TARGET)); 
hf = newfig(title_fig);

set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.4,0.1,0.4,0.6] ) ;
 
Nrow = Nrad; 
Ncol = Ntrial+1;
mergeElement = {};
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 

for mm = 1:Ntrial
Nframe = T{mm}.Nframe; 
id = Nframe; 

h_quiv1 = gobjects(1,Ncirc);
h_quiv2 = gobjects(1,Ncirc);

for iframeMap = 1:length(id)
    ptFC_cur = T{mm}.FC.xyz(:,:,id(iframeMap)); %r.FaceCentroids{id(iframeMap)};
    FC_epc1 = T{mm}.FC.Epc1.cum(:,id(iframeMap)); D1 = T{mm}.FC.Epc1.vec;
    FC_epc2 = T{mm}.FC.Epc2.cum(:,id(iframeMap)); D2 = T{mm}.FC.Epc2.vec;
    for kk =1:Nrad   
        ax_idx = mm+(kk-1)*Ncol;
        mask_rad_current = T{mm}.mask_rad{kk};
        scatter(h_ax(ax_idx),ptFC_cur(:,1),...
            ptFC_cur(:,2),...
            1,'Marker','.',...
            'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
        
        for xx = 1:Ncirc
            mask_circ_current = T{mm}.mask_circ{xx};
            %Fill with NaN outside region
            current_mask = mask_rad_current(:,1)&mask_circ_current(:,1);
            xyz = ptFC_cur(current_mask,:);
            FC_epc1_region = FC_epc1(current_mask,:);
            FC_epc2_region = FC_epc2(current_mask,:);
            D1_region = D1(current_mask,:);
            D2_region = D2(current_mask,:);
            
            % fit ellipse from last contact contour
            [~,centerEl] = myfitellipse(T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2));
            centerEl = [centerEl,zplane_mean];
            pnow = mean(xyz,1,'omitnan'); %centroid

            vecunit_now = pnow-centerEl;
            product_vec1 = dot(repmat(vecunit_now(1:2),length(D1_region),1),D1_region(:,1:2),2);
            product_vec2 = dot(repmat(flip(vecunit_now(1:2)),length(D2_region),1),D2_region(:,1:2),2);
            orientation1 = (product_vec1>0)*2-1; %positive = +1 / negative = -1
            orientation2 = (product_vec2>0)*2-1; %positive = +1 / negative = -1
            D1_region(:,1:2) = D1_region(:,1:2).*orientation1;
            D2_region(:,1:2) = D2_region(:,1:2).*orientation2;
            Dsnow1 = median(D1_region,1,'omitnan').*abs(median(FC_epc1_region,1,'omitnan'))/2; %scaled vector 
            Dsnow2 = median(D2_region,1,'omitnan').*abs(median(FC_epc2_region,1,'omitnan'))/2; %scaled vector                
            
            %plot X-Y points of current region 
            scatter(h_ax(ax_idx),xyz(:,1),xyz(:,2),5,'Marker',marker_set,...
                'MarkerEdgeColor',color_set_cell{xx},'MarkerFaceColor',color_set_cell{xx});
            %plot principales directions projected onto 2D plane 
            h_quiv1(xx) = quiver(h_ax(ax_idx),pnow(1)-Dsnow1(:,1)/2,pnow(2)-Dsnow1(:,2)/2,Dsnow1(:,1),Dsnow1(:,2),0,...
                'Color',color_contraction,'ShowArrowHead','off','AutoScale','off','LineWidth',2);
            h_quiv2(xx) = quiver(h_ax(ax_idx),pnow(1)-Dsnow2(:,1)/2,pnow(2)-Dsnow2(:,2)/2,Dsnow2(:,1),Dsnow2(:,2),0,...
                'Color',color_elongation,'ShowArrowHead','off','AutoScale','off','LineWidth',2);
        end
        % plot contour 
        plot(h_ax(ax_idx),T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2),'k-','LineWidth',1,'Color',[1 1 1].*0.3);
        
        % title(h_ax(ax_idx),region_name_set{xx});
        % xlabel(h_ax(ax_idx),'x'); ylabel(h_ax(ax_idx),'y');
        if mm == 1 
            text(h_ax(ax_idx),min(ptFC_cur(:,1),[],1,'omitnan'),mean(ptFC_cur(:,2),1,'omitnan'),0,T{mm}.name_rad{kk},...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',15,...
                'FontWeight','bold');
        end
        
        if kk == 1
            text(h_ax(ax_idx),mean(ptFC_cur(:,1),1,'omitnan'),min(ptFC_cur(:,2),[],1,'omitnan'),0,sprintf('trial #%d',T{mm}.trial_nbr),...
            'Rotation',0,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15,...
            'FontWeight','bold');
        end
        xlim(h_ax(ax_idx),[-10,10]);
        ylim(h_ax(ax_idx),[-5,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');    
%         set(findall(h_ax(ax_idx), 'type', 'text'), 'visible', 'on')

        % set quiver plot to front of the axis.
        uistack(h_quiv1,'top');
        uistack(h_quiv2,'top');
    end
end
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off','visible','off','Color','None');

% Legend - 1°
ax_idx = 2*Ncol; 
legend_data = cell(1,2);
qw = cell(1,2);
qw{1} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_contraction);
legend_data{1} = sprintf("%s",'Min ~ Contraction');
qw{2} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_elongation);
legend_data{2} = sprintf("%s",'Max ~ Elongation');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','northeast','Interpreter','none');
hleg.Box = 'off';
hleg.NumColumns = 1;
hleg.Title.String = 'Princ. dir.'; 
hleg.Orientation = 'horizontal';
set(hleg,'FontSize',12); 

% Legend - 2°
ax_idx = Ncol;%2*NplotTimeSize; 
legend_data = cell(1,Nrad+1);
qw = cell(1,Nrad+1);
for kk = 1:Nrad
    qw{kk} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_set_cell{kk});
    legend_data{kk} = sprintf("%s",lgd_name{kk});
end
qw{kk+1} = plot(h_ax(ax_idx),NaN,NaN,'-','Color',[1 1 1].*0.4);
legend_data{kk+1} = sprintf("%s",'Last contour');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','best');
hleg.Box = 'off';
hleg.NumColumns = 1;
hleg.Title.String = 'Region'; 
set(hleg,'FontSize',12); 

%% Principal directions
%--------------------------------------------------------------------------
% Vector visualisation 
%--------------------------------------------------------------------------
face_target = 1:2; 
trial_target = 1; 
% settings : colormap, legend name,...
C = colororder(autumn(Ncirc)); %parula
color_set_array = C(1:Ncirc,:);
marker_set = 'o'; 
lgd_name = T{1}.name_circ;
color_set_cell = mat2cell(color_set_array,repelem(1,size(color_set_array,1)),[3]);
color_face = [[0 0 1];[0 1 0]];

% Figure and subplot initialization 
title_fig = sprintf("%s - Principal directions illu = %d",SUBJECT,T{trial_target}.trial_nbr); 
hf = newfig(title_fig);
set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1943    0.4296    0.4182    0.3528] ) ;
 
Nrow = 1; 
Ncol = 3;%Ntrial+3;
mergeElement = {};
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 


p_all = NaN(Ntrial,length(face_target),(Nrad-1)*Ncirc,3); 
Dsunit_all = NaN(Ntrial,length(face_target),(Nrad-1)*Ncirc,3); 
for mm = 1:Ntrial
    Nframe = T{mm}.Nframe;
    id = Nframe;
    
    for iframeMap = 1:length(id)
        ptFC_cur = T{mm}.FC.xyz(:,:,id(iframeMap)); %r.FaceCentroids{id(iframeMap)};
        for iface = face_target
            ax_idx = 1;%+(iface-1)*Ncol;
            face_measure_ii = FACEMEASURE{iface};
            FC_cur = T{mm}.FC.(face_measure_ii).cum(:,id(iframeMap)); 
            D = T{mm}.FC.(face_measure_ii).vec;
            % fit ellipse from last contact contour
            [CAfit,centerEl] = myfitellipse(T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2));
            %
            CAfit_1 = (CAfit - centerEl)*(1/3)+centerEl;
            CAfit_2 = (CAfit - centerEl)*(2/3)+centerEl;
            CAfit_3 = CAfit;
            
            zplane_mean = mean(T{mm}.ca{end}(:,3)); 
            CAfit_1 = [CAfit_1,repelem(zplane_mean,size(CAfit,1),1)];
            CAfit_2 = [CAfit_2,repelem(zplane_mean,size(CAfit,1),1)];
            CAfit_3 = [CAfit_3,repelem(zplane_mean,size(CAfit,1),1)];
            centerEl = [centerEl,zplane_mean];
            
            ll = 0; 
            for kk =1:Nrad-1
                mask_rad_current = T{mm}.mask_rad{kk};
                h_quiv = gobjects(1,Ncirc);
                h_quivpt = gobjects(1,Ncirc);
                for xx = 1:Ncirc
                    ll = ll+1; 
                    mask_circ_current = T{mm}.mask_circ{xx};
                    %Fill with NaN outside region
                    current_mask = mask_rad_current(:,1)&mask_circ_current(:,1);
                    xyz = ptFC_cur(current_mask,:);
                    FC_region = FC_cur(current_mask,:);
                    vec_region = D(current_mask,:);
                    
                    pnow = mean(xyz,1,'omitnan'); %centroid
                    
                    vecunit_now = pnow-centerEl;
                    if iface == 1
                        product_vec = dot(repmat(vecunit_now(1:2),length(vec_region),1),vec_region(:,1:2),2);
                    elseif iface == 2
                        product_vec = dot(repmat(flip(vecunit_now(1:2)),length(vec_region),1),vec_region(:,1:2),2);
                    end
                    orientation = (product_vec>0)*2-1; %positive = +1 / negative = -1
                    vec_region(:,1:2) = vec_region(:,1:2).*orientation;
                    
                    Dsnow = median(vec_region,1,'omitnan'); %.*median(FC_region,1,'omitnan'); %scaled vector                 
                    Dsunitnow = Dsnow/vecnorm(Dsnow,2,2)*2; %*median(FC_region,1,'omitnan'); 
                    
                    p_all(mm,iface,ll,:) = pnow';
                    Dsunit_all(mm,iface,ll,:) = Dsunitnow';
                    if mm == trial_target
                    %plot X-Y points of current region 
                    k = boundary(xyz(:,1),xyz(:,2),0.3); 
                    plot(h_ax(ax_idx),xyz(k,1),xyz(k,2),'k-');
                    %plot principales directions projected onto 2D plane
                    h_quiv(xx) = quiver(h_ax(ax_idx),pnow(1)-Dsunitnow(:,1)/2,pnow(2)-Dsunitnow(:,2)/2,Dsunitnow(:,1),Dsunitnow(:,2),0,...
                        'Color',color_face(iface,:),'ShowArrowHead','off','AutoScale','off','LineWidth',3);
                    h_quivpt(xx) = plot(h_ax(ax_idx),pnow(1),pnow(2),'s','Color',[0.8500 0.3250 0.0980],...
                        'MarkerSize',3,'MarkerFaceColor',[0.8500 0.3250 0.0980]); 
                    end
                end
            end
            
            if mm == trial_target
            h_mapscatter = scatter(h_ax(ax_idx),ptFC_cur(:,1),...
                ptFC_cur(:,2),...
                2,'Marker','.',...
                'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
            % plot contour area bin
            plot(h_ax(ax_idx),centerEl(:,1),centerEl(:,2),'k.','MarkerSize',20);
            if iface == 1
            text(h_ax(ax_idx),min(ptFC_cur(:,1),[],1,'omitnan')-3,mean(ptFC_cur(:,2),1,'omitnan'),0,...
                sprintf('Principal\ndirections'),...
                'Rotation',90,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','center',...
                'FontSize',18,...
                'FontWeight','bold');
            end
            xlim(h_ax(ax_idx),[-12,12]);
            ylim(h_ax(ax_idx),[-5,15]);
            set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');
            %set quiver plot to front of the axis.
            uistack(h_quiv,'top');
            uistack(h_quivpt,'top');
            uistack(h_mapscatter,'bottom');
            end
            
%             jump = 1; 
%             quiver(h_ax(ax_idx),ptFC_cur(1:jump:end,1)-D(1:jump:end,1)/2,ptFC_cur(1:jump:end,2)-D(1:jump:end,2)/2,D(1:jump:end,1)/10,D(1:jump:end,2)/2,0,...
%                         'Color',color_face(iface,:),'ShowArrowHead','off','AutoScale','off','LineWidth',1);
        end
    end
end
orientation_face = ["Meridional","Circumferencial"];
for iface = face_target
    ax_idx = Ncol-2+iface;
    if Ntrial~=1
        pnow = squeeze(mean(squeeze(p_all(:,iface,:,:)),1,'omitnan'));
        Dsunitnow = squeeze(mean(squeeze(Dsunit_all(:,iface,:,:)),1,'omitnan'));
    else
        pnow = squeeze(p_all(:,iface,:,:));
        Dsunitnow = squeeze(Dsunit_all(:,iface,:,:));
    end
    vecunit_now = (pnow(end-4:end,:)-centerEl)*1.2;
        
            h_mapscatter = scatter(h_ax(ax_idx),ptFC_cur(:,1),...
                ptFC_cur(:,2),...
                3,'Marker','.',...
                'MarkerEdgeColor',[1 1 1].*0.2,'MarkerFaceColor',[1 1 1].*0.2);

    %plot principales directions projected onto 2D plane
    h_quiv = quiver(h_ax(ax_idx),pnow(:,1)-Dsunitnow(:,1)/2,pnow(:,2)-Dsunitnow(:,2)/2,Dsunitnow(:,1),Dsunitnow(:,2),0,...
        'Color',color_face(iface,:),'ShowArrowHead','off','AutoScale','off','LineWidth',3);
    h_quivpt = plot(h_ax(ax_idx),pnow(:,1),pnow(:,2),'s','Color',[0.8500 0.3250 0.0980],...
        'MarkerSize',3,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
    plot(h_ax(ax_idx),centerEl(:,1),centerEl(:,2),'k.','MarkerSize',20);
    if iface == 2
    text(h_ax(ax_idx),mean(ptFC_cur(:,1),1,'omitnan')-15,min(ptFC_cur(:,2),[],1,'omitnan'),0,...
        sprintf('Trial-averaged'),...
        'Rotation',0,...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','center',...
        'FontSize',15,...
        'FontWeight','bold');
    end
    text(h_ax(ax_idx),centerEl(1),centerEl(2)+5,0,...
        sprintf('%s',orientation_face(iface)),...
        'Rotation',0,...
        'VerticalAlignment','top',...
        'HorizontalAlignment','center',...
        'FontSize',15,...
        'FontWeight','bold');

    xlim(h_ax(ax_idx),[-12,12]);
    ylim(h_ax(ax_idx),[-5,15]);
    set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');
    % set quiver plot to front of the axis.
    uistack(h_quiv,'top');
    uistack(h_quivpt,'top');
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off','visible','off','Color','None');

% Legend - 1°
ax_idx = Ncol; 
legend_data = cell(1,2);
qw = cell(1,2);
qw{1} = fill(h_ax(ax_idx),NaN,NaN,color_face(1,:),'EdgeColor','none');
legend_data{1} = sprintf("%s",'Epc1');
qw{2} = fill(h_ax(ax_idx),NaN,NaN,color_face(2,:),'EdgeColor','none');
legend_data{2} = sprintf("%s",'Epc2');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','northeast','Interpreter','none');
hleg.Box = 'off';
hleg.NumColumns = 2;
% hleg.Title.String = 'Princ. dir.'; 
hleg.Orientation = 'horizontal';
set(hleg,'FontSize',12,'Position',[0.5686    0.2120    0.2042    0.0551]); 



exportgraphics(hf,sprintf('fig/%s.%s',title_fig,'png')); 

%--------------------------------------------------------------------------
%% SUMMARY SEGMENTATION
%--------------------------------------------------------------------------

id = [1,91]; 

% settings : colormap, legend name,...
C = colororder(autumn(Ncirc)); %parula
color_set_array = C(1:Ncirc,:);
marker_set = 'o'; 
lgd_name = T{1}.name_circ;
color_set_cell = mat2cell(color_set_array,repelem(1,size(color_set_array,1)),[3]);

% Figure and subplot initialization 
title_fig = sprintf("%s : Map point cloud segmentation summary\n%s",SUBJECT,mat2str(TRIAL_TARGET)); 
hf = newfig(title_fig);

set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.4,0.1,0.4,0.6] ) ;
 
Nrow = 4; 
Ncol = Ntrial+1;
mergeElement = {};
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 

for mm = 1:Ntrial

h_quiv1 = gobjects(1,Ncirc);
h_quiv2 = gobjects(1,Ncirc);

for iframeMap = 1:length(id)
    ptFC_cur = T{mm}.FC.xyz(:,:,id(iframeMap)); %r.FaceCentroids{id(iframeMap)};
    for kk =1:Nrad   
        ax_idx = mm+(iframeMap-1)*Ncol;
        mask_rad_current = T{mm}.mask_rad{kk};
        if kk == 1
            scatter3(h_ax(ax_idx),ptFC_cur(:,1),...
                ptFC_cur(:,2),...
                ptFC_cur(:,3),...
                1,'Marker','.',...
                'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
        end
        for xx = 1:Ncirc
            mask_circ_current = T{mm}.mask_circ{xx};
            %Fill with NaN outside region
            current_mask = mask_circ_current(:,1)&mask_rad_current(:,1);
            xyz = ptFC_cur(current_mask,:);
            
            %plot X-Y points of current region 
            scatter3(h_ax(ax_idx),xyz(:,1),xyz(:,2),xyz(:,3),5,'Marker',marker_set,...
                'MarkerEdgeColor',color_set_cell{xx},'MarkerFaceColor',color_set_cell{xx});
        end
        % plot contour 
        plot3(h_ax(ax_idx),T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2),T{mm}.ca{end}(:,3),'k-','LineWidth',1,'Color',[1 1 1].*0.3);
        
        if kk == 1
            text(h_ax(ax_idx),mean(ptFC_cur(:,1),1,'omitnan'),min(ptFC_cur(:,2),[],1,'omitnan'),0,sprintf('trial #%d',T{mm}.trial_nbr),...
            'Rotation',0,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15,...
            'FontWeight','bold');
        end
        xlim(h_ax(ax_idx),[-10,10]);
        ylim(h_ax(ax_idx),[-5,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');    
%         set(findall(h_ax(ax_idx), 'type', 'text'), 'visible', 'on')

    end
end
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off','visible','off','Color','None');

% Legend - 1°
ax_idx = Ncol;%2*NplotTimeSize; 
legend_data = cell(1,Nrad+1);
qw = cell(1,Nrad+1);
for kk = 1:Ncirc
    qw{kk} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_set_cell{kk});
    legend_data{kk} = sprintf("%s",lgd_name{kk});
end
qw{kk+1} = plot(h_ax(ax_idx),NaN,NaN,'-','Color',[1 1 1].*0.4);
legend_data{kk+1} = sprintf("%s",'Last contour');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','best');
hleg.Box = 'off';
hleg.NumColumns = 1;
hleg.Title.String = 'Region'; 
set(hleg,'FontSize',12); 

%--------------------------------------------------------------------------
% Segmentation RADIAL 
%--------------------------------------------------------------------------

% settings : colormap, legend name,...
C=[0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
color_set_array = C(1:Nrad,:);
marker_set = 'o'; 
lgd_name = T{1}.name_rad;
color_set_cell = mat2cell(color_set_array,repelem(1,size(color_set_array,1)),[3]);

for mm = 1:Ntrial

for iframeMap = 1:length(id)
    ptFC_cur = T{mm}.FC.xyz(:,:,id(iframeMap)); %r.FaceCentroids{id(iframeMap)};
    for xx = 1:Ncirc
        ax_idx = mm+(iframeMap-1)*Ncol+length(id)*Ncol;
        mask_circ_current = T{mm}.mask_circ{xx};
        scatter3(h_ax(ax_idx),ptFC_cur(:,1),...
            ptFC_cur(:,2),...
            ptFC_cur(:,3),...
            1,'Marker','.',...
            'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
        
        
        for kk =1:Nrad-1
            mask_rad_current = T{mm}.mask_rad{kk};
            %Fill with NaN outside region
            current_mask = mask_rad_current(:,1)&mask_circ_current(:,1);
            xyz = ptFC_cur(current_mask,:);
            
            %plot X-Y points of current region 
            scatter3(h_ax(ax_idx),xyz(:,1),xyz(:,2),xyz(:,3),5,'Marker',marker_set,...
                'MarkerEdgeColor',color_set_cell{kk},'MarkerFaceColor',color_set_cell{kk});
        end
        % plot contour 
        plot3(h_ax(ax_idx),T{mm}.ca{end}(:,1),T{mm}.ca{end}(:,2),T{mm}.ca{end}(:,3),'k-','LineWidth',1,'Color',[1 1 1].*0.3);
        
        if xx == 1
            text(h_ax(ax_idx),mean(ptFC_cur(:,1),1,'omitnan'),min(ptFC_cur(:,2),[],1,'omitnan'),0,sprintf('trial #%d',T{mm}.trial_nbr),...
            'Rotation',0,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15,...
            'FontWeight','bold');
        end
        xlim(h_ax(ax_idx),[-10,10]);
        ylim(h_ax(ax_idx),[-5,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90.0000],'visible','off');    

    end
end
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off','visible','off','Color','None');

% Legend - 2°
ax_idx = 2*Ncol;%2*NplotTimeSize; 
legend_data = cell(1,Nrad+1);
qw = cell(1,Nrad+1);
for kk = 1:Nrad
    qw{kk} = fill(h_ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],color_set_cell{kk});
    legend_data{kk} = sprintf("%s",lgd_name{kk});
end
qw{kk+1} = plot(h_ax(ax_idx),NaN,NaN,'-','Color',[1 1 1].*0.4);
legend_data{kk+1} = sprintf("%s",'Last contour');
hleg = legend(h_ax(ax_idx),[qw{:}],legend_data,'Location','best','Interpreter','none');
hleg.Box = 'off';
hleg.NumColumns = 1;
hleg.Title.String = 'Region'; 
set(hleg,'FontSize',12); 


%--------------------------------------------------------------------------
%% SUMMARY CONTACT BORDER
%--------------------------------------------------------------------------
color_annulus = 'b';
marker_set = '.';
% Figure and subplot initialization 
title_fig = sprintf("%s : Map point cloud contact border summary\n%s",SUBJECT,mat2str(TRIAL_TARGET)); 
hf = newfig(title_fig);

set(hf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.4,0.1,0.4,0.6] ) ;

id = [1,15,20,21,22,23,24,25,Nframe]; 
Nrow = Ntrial; 
Ncol = length(id)+1;
mergeElement = {};
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 

for iframeMap = 1:length(id)
for mm = 1:Ntrial
    ptFC_cur = T{mm}.FC.xyz(:,:,id(iframeMap)); %r.FaceCentroids{id(iframeMap)};
        ax_idx = iframeMap+(mm-1)*Ncol;
        scatter3(h_ax(ax_idx),ptFC_cur(:,1),...
            ptFC_cur(:,2),...
            ptFC_cur(:,3),...
            1,'Marker','.',...
            'MarkerEdgeColor',[1 1 1].*0.5,'MarkerFaceColor',[1 1 1].*0.5);
        
        mask_in_annulus_current = T{mm}.mask_in_annulus(:,id(iframeMap)); %&~T{mm}.mask_in(:,id(iframeMap)); 
        mask_in_border_current = T{mm}.mask_in_border(:,id(iframeMap));
        
        %Fill with NaN outside region
        current_mask = mask_in_annulus_current;
        xyz = ptFC_cur(current_mask,:);
        
        ca_now = T{mm}.ca{id(iframeMap)};
        if length(ca_now)>3
            ca_now_fitted = myfitellipse(ca_now(:,1),ca_now(:,2)); 
            ca_now_fitted = [ca_now_fitted,repelem(mean(ca_now(:,3),'omitnan'),length(ca_now_fitted),1)];
        else
            ca_now_fitted = [NaN,NaN,NaN]; 
        end
        
        %plot X-Y points of current region 
        scatter3(h_ax(ax_idx),xyz(:,1),xyz(:,2),xyz(:,3),2,'Marker',marker_set,...
            'MarkerEdgeColor',color_annulus,'MarkerFaceColor',color_annulus);
%         % Point on the boundary 
%         scatter3(h_ax(ax_idx),ptFC_cur(current_mask,1),ptFC_cur(current_mask,2),ptFC_cur(current_mask,3),15,'r');
        % plot contour 
        plot3(h_ax(ax_idx),ca_now(:,1),ca_now(:,2),ca_now(:,3),'k-','LineWidth',2);
        % plot fitted ellipse 
        plot3(h_ax(ax_idx),ca_now_fitted(:,1),ca_now_fitted(:,2),ca_now_fitted(:,3),'-','Color',[0.8500 0.3250 0.0980],'LineWidth',3);
        
        if mm == 1
            text(h_ax(ax_idx),mean(ptFC_cur(:,1),1,'omitnan'),min(ptFC_cur(:,2),[],1,'omitnan'),0,sprintf('frame #%d',id(iframeMap)),...
            'Rotation',0,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15,...
            'FontWeight','bold');
        end
        
        if iframeMap == 1
            text(h_ax(ax_idx),min(ptFC_cur(:,1),[],1,'omitnan'),mean(ptFC_cur(:,2),1,'omitnan'),0,sprintf('trial #%d',T{mm}.trial_nbr),...
            'Rotation',90,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15,...
            'FontWeight','bold');
        end
        xlim(h_ax(ax_idx),[-10,10]);
        ylim(h_ax(ax_idx),[-5,15]);
        set(h_ax(ax_idx),'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90],'visible','off');    
%         set(findall(h_ax(ax_idx), 'type', 'text'), 'visible', 'on')

end
end
% Set global axis parameters
set(h_ax,'FontSize',12,'Clipping','off','visible','off','Color','None');


%% COMPUTATION annulus - time dependent variables 
fprintf('Computation annulus - time varying variables...'); 
Nframe_max_all = 0; 
for mm = 1:Ntrial
    Nframe_cur = T{mm}.Nframe; 
    if Nframe_cur > Nframe_max_all 
        Nframe_max_all = Nframe_cur; 
    end
end
FC_rate_val_in_annulus_all = NaN(Nframe_max_all,Nface,Ntrial); 
FC_cum_val_in_annulus_all = NaN(Nframe_max_all,Nface,Ntrial); 
FC_rate_val_out_annulus_all = NaN(Nframe_max_all,Nface,Ntrial); 
FC_cum_val_out_annulus_all = NaN(Nframe_max_all,Nface,Ntrial); 
FC_rate_val_annulus_all = NaN(Nframe_max_all,Nface,Ntrial); 
FC_cum_val_annulus_all = NaN(Nframe_max_all,Nface,Ntrial); 
for mm = 1:Ntrial
    Nframe = T{mm}.Nframe; 
    Ndata = length(T{mm}.FC.(FACEMEASURE{1}).cum); 
    mask_in_annulus_current = T{mm}.mask_in_annulus; 
    mask_in_current = T{mm}.mask_in(:,id(iframeMap));
    
    for iface = 1:2;%Nface
        faceMeasure_ii = FACEMEASURE{iface}; 
        FC_cum = T{mm}.FC.(faceMeasure_ii).cum; 
        FC_rate = T{mm}.FC.(faceMeasure_ii).rate; 
        
        %current var
        FC_cum_current = FC_cum;
        FC_rate_current = FC_rate; 
        idx_frontal_target = T{mm}.mask_circ{1}(:,1)|T{mm}.mask_circ{2}(:,1)|T{mm}.mask_circ{3}(:,1)|T{mm}.mask_circ{4}(:,1)|T{mm}.mask_circ{5}(:,1); 
        current_mask = mask_in_annulus_current&idx_frontal_target; 
        FC_cum_current(~current_mask) = NaN;
        FC_rate_current(~current_mask) = NaN;
        
        %var computation
        if strcmp(faceMeasure_ii,'Area')
            FC_rate_val_current = sum(FC_cum_current,1,'omitnan');
            FC_cum_val_current = sum(FC_cum_current,1,'omitnan');
        else
            FC_rate_val_current = median(FC_rate_current,1,'omitnan');
            FC_cum_val_current = median(FC_cum_current,1,'omitnan');
        end
        FC_rate_val_annulus_all(1:Nframe,iface,mm) = FC_rate_val_current; 
        FC_cum_val_annulus_all(1:Nframe,iface,mm) = FC_cum_val_current; 
        
        %current var
        FC_cum_current = FC_cum;
        FC_rate_current = FC_rate; 
        current_mask = mask_in_annulus_current&mask_in_current&idx_frontal_target; 
        FC_cum_current(~current_mask) = NaN;
        FC_rate_current(~current_mask) = NaN;
        
        %var computation
        if strcmp(faceMeasure_ii,'Area')
            FC_rate_val_current = sum(FC_cum_current,1,'omitnan');
            FC_cum_val_current = sum(FC_cum_current,1,'omitnan');
        else
            FC_rate_val_current = median(FC_rate_current,1,'omitnan');
            FC_cum_val_current = median(FC_cum_current,1,'omitnan');
        end
        FC_rate_val_in_annulus_all(:,iface,mm) = FC_rate_val_current; 
        FC_cum_val_in_annulus_all(:,iface,mm) = FC_cum_val_current; 
        
        %current var
        FC_cum_current = FC_cum;
        FC_rate_current = FC_rate; 
        idx_frontal_target = T{mm}.mask_circ{1}(:,1)|T{mm}.mask_circ{2}(:,1)|T{mm}.mask_circ{3}(:,1)|T{mm}.mask_circ{4}(:,1)|T{mm}.mask_circ{5}(:,1); 
        current_mask = mask_in_annulus_current&~mask_in_current&idx_frontal_target; 
        FC_cum_current(~current_mask) = NaN;
        FC_rate_current(~current_mask) = NaN;
        
        %var computation
        if strcmp(faceMeasure_ii,'Area')
            FC_rate_val_current = sum(FC_cum_current,1,'omitnan');
            FC_cum_val_current = sum(FC_cum_current,1,'omitnan');
        else
            FC_rate_val_current = median(FC_rate_current,1,'omitnan');
            FC_cum_val_current = median(FC_cum_current,1,'omitnan');
        end
        FC_rate_val_out_annulus_all(1:Nframe,iface,mm) = FC_rate_val_current; 
        FC_cum_val_out_annulus_all(1:Nframe,iface,mm) = FC_cum_val_current; 
        
    end
    %output var
    T{mm}.FC_rate_val_annulus = FC_rate_val_annulus_all(:,:,mm);
    T{mm}.FC_cum_val_annulus = FC_cum_val_annulus_all(:,:,mm);
    T{mm}.FC_rate_val_in_annulus = FC_rate_val_in_annulus_all(:,:,mm);
    T{mm}.FC_cum_val_in_annulus = FC_cum_val_in_annulus_all(:,:,mm);
    T{mm}.FC_rate_val_out_annulus = FC_rate_val_out_annulus_all(:,:,mm);
    T{mm}.FC_cum_val_out_annulus = FC_cum_val_out_annulus_all(:,:,mm);
    
end

fprintf('done\n'); 

%% POLE TO EDGES ANNULUS

fprintf(' Multiple segmentation logic...')

Nrad_out_set = 7; 
Nrad_in_set = 10; 
rad_combined_number = 3; 
Nrad_segmentation = Nrad_in_set+Nrad_out_set; 
% FACEMEASURE{4}='curv'; 
Nface = length(FACEMEASURE); 
trial_target = 1:Ntrial; 
p_in_threshold_inside_logic = 0.10; 
% frame_target = 95;%0.5N : 21 % 1N : 30 % 2N : 43 % 5N : 76
list_force_target = [0.00,0.01,0.02,0.03,0.04,0.05,0.10,0.20,0.50,1:1:5];
% settings : colormap, legend name,...
C=[0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840;...
    0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840;...
    0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840;...
    0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
% C = [1 1 1].*1-linspace(0.1,1,Nrad_segmentation)';
% C = colormap(gray(round(Nrad_segmentation*(1.2)))); 
color_set_array = C(1:Nrad_segmentation,:);
idx_IN_rad_all = cell(1,length(trial_target));
FCmat_rate = cell(1,Ntrial);
FCmat_rate_last_all = cell(1,Ntrial);
FCmat_cum = cell(1,Ntrial);
FCmat_xyz = cell(1,Ntrial);
FCmat_xyz_all = cell(1,Ntrial); 
FCmat_cum_last_all = cell(1,Ntrial);
npoints_compute_all = cell(1,Ntrial); 
logic_inside_rad_all = cell(1,Ntrial);
Area_all = cell(1,Ntrial); 

for mm = trial_target
    Nrad_out = Nrad_out_set; 
    Nrad_in = Nrad_in_set; 
    Nrad_segmentation = Nrad_in+Nrad_out; 
    
%     Nframe = T{mm}.Nframe; 
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>5,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    frame_target = list_frame_num(end); 
    
    PCtarget = T{mm}.FC.xyz(:,:,frame_target);
    CAtarget = T{mm}.ca{frame_target};
    idx_frontal_target = T{mm}.mask_circ{1}(:,1)|T{mm}.mask_circ{2}(:,1)|T{mm}.mask_circ{3}(:,1)|T{mm}.mask_circ{4}(:,1)|T{mm}.mask_circ{5}(:,1); 

    idx_rad_target = T{mm}.mask_circ{2}(:,1)|T{mm}.mask_circ{3}(:,1)|T{mm}.mask_circ{4}(:,1); 
    
    % fit ellipse from last contact contour
    [CAfit,centerEl,Rx,Ry] = myfitellipse(CAtarget(:,1),CAtarget(:,2));
    z_height = mean(CAtarget(:,3),'omitnan');
    %
    
    %define ellipses inside : using a fraction of the last contour
    CAfit_rad = cell(1,Nrad_segmentation);
    idx_IN_rad_raw = cell(1,Nrad_segmentation);
    idx_IN_allprevious = false(length(PCtarget),1);
    for irad = 1:Nrad_in    
        if irad ~= 1
            idx_IN_allprevious = idx_IN_allprevious|idx_IN_rad_raw{irad-1};
        end
        CAfit_rad{irad} = (CAfit - centerEl)*(irad/Nrad_in)+centerEl;    
        idx_IN_current = inpolygon(PCtarget(:,1),PCtarget(:,2),CAfit_rad{irad}(:,1),CAfit_rad{irad}(:,2));
        idx_IN_rad_raw{irad} = idx_IN_current & ~idx_IN_allprevious & idx_rad_target;
    end
    Nrad_border = 1;
    for irad = Nrad_in+1:Nrad_in+Nrad_border
        idx_IN_allprevious = idx_IN_allprevious|idx_IN_rad_raw{irad-1};
        CAfit_rad{irad} = (CAfit - centerEl)*(irad/Nrad_in)+centerEl;    
        idx_IN_current = inpolygon(PCtarget(:,1),PCtarget(:,2),CAfit_rad{irad}(:,1),CAfit_rad{irad}(:,2));
        idx_IN_rad_raw{irad} = idx_IN_current & ~idx_IN_allprevious & idx_rad_target;
    end
    irad_last = irad; 
    %define ellipses outside : 
    mask_in_target = T{mm}.mask_in(:,frame_target);
    p_origin = mean(PCtarget(mask_in_target,:),'omitnan');
    z_step  = 0.5*Rx/Nrad_in; % increment of distance in mm
    ii = 0; 
    for irad = irad_last+1:Nrad_segmentation
        ii = ii + 1; 
        idx_IN_current = PCtarget(:,3)<(p_origin(3)+ii*z_step);
        idx_IN_allprevious = idx_IN_allprevious|idx_IN_rad_raw{irad-1};
        idx_IN_rad_raw{irad} = idx_IN_current & ~idx_IN_allprevious & idx_rad_target;
    end
        
    %Merge (#margin) number of radial segmentation at the origin without
    %the very center where the orientation of deformation is not well
    %defined.
    idx_IN_rad_first = false(length(PCtarget),1); 
    for irad = 2:rad_combined_number
        idx_IN_rad_first = idx_IN_rad_first|idx_IN_rad_raw{irad};
    end
    idx_IN_rad = cell(1,Nrad_segmentation-rad_combined_number+1);
    idx_IN_rad{1} = idx_IN_rad_first;
    idx_IN_rad(2:end) = idx_IN_rad_raw(1+rad_combined_number:end);
    
    Nrad_segmentation = Nrad_segmentation-rad_combined_number+1;
    Nrad_in = Nrad_in-rad_combined_number+1; 
    
    CAfit =[CAfit,repelem(mean(CAtarget(:,3),'omitnan'),length(CAfit),1)];
    PCcurrent = T{mm}.FC.xyz(:,:,frame_target);

    % median val
    Ndata = length(PCcurrent); 
    Nface = length(FACEMEASURE); 
    FCmat_cum_last_all_cur = NaN(Ndata,Nface*(Nrad_segmentation)); 
    FCmat_cum_cur = NaN(Nframe,Nface*(Nrad_segmentation)); 
    FCmat_rate_cur = NaN(Nframe,Nface*(Nrad_segmentation)); 
    FCmat_xyz_cur = NaN(Nframe,3,(Nrad_segmentation)); 
    FCmat_xyz_all_cur = NaN(Ndata,Nframe,3,(Nrad_segmentation)); 
    npoints_compute_cur = NaN(Nrad_segmentation,1); 
    logic_inside_rad_cur = false(Nframe,Nrad_segmentation);
    iface = 0; 
    for face_num = 1:Nface
        iface = iface + 1; 
        faceMeasure_ii = FACEMEASURE{face_num};
        FC_cum = T{mm}.FC.(faceMeasure_ii).cum; 
        FC_rate = T{mm}.FC.(faceMeasure_ii).rate; 
        FC_xyz = T{mm}.FC.xyz;
        for irad = 1:Nrad_segmentation
            %current var
            FC_cum_current = FC_cum;
            FC_rate_current = FC_rate;
            FC_xyz_current = FC_xyz;
            mask_rad_current = idx_IN_rad{irad};
            current_mask = mask_rad_current;
            %
            FC_cum_current(~current_mask,:) = NaN;
            FC_rate_current(~current_mask,:) = NaN;
            FC_xyz_current(~current_mask,:,:) = NaN;
            %
            if iface == 1
                FCmat_xyz_all_cur(:,:,:,irad) = permute(FC_xyz_current,[1 3 2]);  
                FCmat_xyz_cur(:,:,irad) = squeeze(median(FC_xyz_current,1,'omitnan'))';       
            end
%             if iface == 5 || iface == 6
%                 FCmat_cum_cur(:,irad+(5)*Nrad_segmentation) = FCmat_cum_cur(:,irad+(5)*Nrad_segmentation)+median(FC_cum_current,1,'omitnan')'; 
%                 FCmat_rate_cur(:,irad+(5)*Nrad_segmentation) = FCmat_rate_cur(:,irad+(5)*Nrad_segmentation)+median(FC_rate_current,1,'omitnan')'; 
%                 FCmat_cum_last_all_cur(:,irad+5*Nrad_segmentation) = FCmat_cum_last_all_cur(:,irad+5*Nrad_segmentation)+FC_cum_current(:,frame_target); %all data at frame_target
%             end
            npoints_compute_cur(irad) = sum(current_mask);
            FCmat_cum_cur(:,irad+(iface-1)*Nrad_segmentation) = median(FC_cum_current,1,'omitnan'); 
            FCmat_rate_cur(:,irad+(iface-1)*Nrad_segmentation) = median(FC_rate_current,1,'omitnan'); 
            FCmat_cum_last_all_cur(:,irad+(iface-1)*Nrad_segmentation) = FC_cum_current(:,frame_target); %all data at frame_target
            %
        end
    end
    
    %inside logic 
    FC_xyz = T{mm}.FC.xyz;
    mask_in = T{mm}.mask_in;
    A = zeros(Nframe,1); 
    for iframe = 1:Nframe
        FC_xyz_current = FC_xyz(:,:,iframe);     
        CA_current = T{mm}.ca{iframe};

        if ~isempty(CA_current) && length(CA_current) > 3
            ca_border_fitted = myfitellipse(CA_current(:,1),CA_current(:,2));
        else
            ca_border_fitted = [NaN,NaN];
        end
        
        for irad = 1:Nrad_segmentation
            FC_xyz_rad_current = FC_xyz_current;
            mask_rad_current = idx_IN_rad{irad};
            current_mask = mask_rad_current;
            FC_xyz_rad_current(~current_mask,:) = NaN;
            mask_inside = inpolygon(FC_xyz_rad_current(:,1),FC_xyz_rad_current(:,2),ca_border_fitted(:,1),ca_border_fitted(:,2));
            logic_inside_rad_cur(iframe,irad) = sum(mask_inside)>p_in_threshold_inside_logic*sum(current_mask);
        end
        
        %Area computation
        if ~any(isnan(ca_border_fitted))
            [x,y] = poly2cw(ca_border_fitted(:,1),ca_border_fitted(:,2));
            A(iframe) = polyarea(x,y);
            if abs(A(iframe))>500
                A(iframe) = NaN; 
            end
        end
    end
    
    %output variables
    FCmat_cum{mm} = FCmat_cum_cur;
    FCmat_rate{mm} = FCmat_rate_cur;
    FCmat_xyz{mm} = FCmat_xyz_cur;
    FCmat_xyz_all{mm} = FCmat_xyz_all_cur;
    FCmat_cum_last_all{mm} = FCmat_cum_last_all_cur;
    npoints_compute_all{mm} = npoints_compute_cur; 
    logic_inside_rad_all{mm} = logic_inside_rad_cur;
    idx_IN_rad_all{mm} = idx_IN_rad;
    Area_all{mm} = A; 
end
clear FCmat_cum_cur FCmat_rate_cur FCmat_xyz_cur FCmat_xyz_all_cur FCmat_cum_last_all_cur logic_inside_rad_cur idx_IN_rad
fprintf('done\n');


%% save output in a table and in a structure 
list_force_target = [0.0,0.025,0.05,0.10,0.15,0.20,0.30,0.40,0.5:0.25:NFCOND_SET]; 
file_table_target = fullfile("analysis",sprintf('%s_table_normal_loading.mat',SUBJECT));
trial_target_list = 1:Ntrial; 
%file loading 
S = dir(file_table_target);
if isempty(S)
    disp('file not found, new table will be created');
    table_normal_loading = table(); 
else 
    disp('file found, table loaded');
    table_normal_loading = load(fullfile(S.folder,S.name)); table_normal_loading = table_normal_loading.table_normal_loading;
end
for mm = trial_target_list
    Nframe = T{mm}.Nframe; 
    spddxlcond_cur = trial_info.Spddxlcond(T{mm}.trial_nbr == trial_info.TrialNbr);
    nfcond_cur = trial_info.Nfcond(T{mm}.trial_nbr == trial_info.TrialNbr);
    
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>nfcond_cur*0.99,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    
    Subject = repelem(SUBJECT,length(list_frame_num),1); 
    Phase = repelem(PHASE,length(list_frame_num),1); 
    Material =  repelem(FRcond{MATERIAL_NBR},length(list_frame_num),1);
    TrialNbr =  repelem(T{mm}.trial_nbr,length(list_frame_num),1);
    NFcond =  repelem(nfcond_cur,length(list_frame_num),1); 
    SPDcond =  repelem(spddxlcond_cur,length(list_frame_num),1);
    FrNbr = list_frame_num'; 
    NF_raw = T{mm}.nf(list_frame_num); 
    NF_formatted = list_force_target'; 
    TF_raw = T{mm}.tf(list_frame_num); 
    Area = Area_all{mm}(list_frame_num);
    
    %-------------
    %to be added 
    angle_cur = T{mm}.angle; 
    Attack_angle = repelem(angle_cur,length(list_frame_num),1);
    %-------------
    
    Ttemp = table(); 
    for iseg = 1:Nrad_segmentation
        SegFactor = repelem(round(iseg/Nrad_in*100)/100,length(list_frame_num),1); 
        
        logic_in = logic_inside_rad_all{mm}(list_frame_num,iseg); 
        
        %-------------
        %to be added 
        npc_computed_cur = npoints_compute_all{mm}(iseg);
        Npc_computed = repelem(npc_computed_cur,length(list_frame_num),1);
        angle_cur = T{mm}.angle; 
        %-------------
        if npc_computed_cur>5
            iface = 1; 
            Epc1_cum = FCmat_cum{mm}(list_frame_num,iseg+(iface-1)*Nrad_segmentation);
            Epc1_rate = FCmat_rate{mm}(list_frame_num,iseg+(iface-1)*Nrad_segmentation);
            iface = 2; 
            Epc2_cum = FCmat_cum{mm}(list_frame_num,iseg+(iface-1)*Nrad_segmentation);
            Epc2_rate = FCmat_rate{mm}(list_frame_num,iseg+(iface-1)*Nrad_segmentation);
            iface = 3; 
            DispZ = FCmat_cum{mm}(list_frame_num,iseg+(iface-1)*Nrad_segmentation);
            iface = 4; 
            Curv = FCmat_cum{mm}(list_frame_num,iseg+(iface-1)*Nrad_segmentation);
            Ttemp_low = table(Subject,Phase,Material,TrialNbr,NFcond,SPDcond,FrNbr,SegFactor,NF_raw,NF_formatted,TF_raw,Area,...
                Epc1_cum,Epc1_rate,Epc2_cum,Epc2_rate,DispZ,Curv,logic_in); 
            Ttemp = [Ttemp;Ttemp_low];  
        end
    end
    
    %Replace table by the temporary one if it is the first iteration
    %or else :
    %    if trial is already inside the table,
    %       replace the old values by newly processed one.
    %    else,
    %       concatenate the table
    %
    if isempty(table_normal_loading)
        table_normal_loading = Ttemp;
    else
        idx_trial = table_normal_loading.TrialNbr==unique(Ttemp.TrialNbr);
        if size(table_normal_loading,2) ~= size(Ttemp,2)
            fprintf('Tables have different number of variables\n'); 
            error('Suppress the last table !'); 
        end
        if any(ismember(table_normal_loading.TrialNbr,Ttemp.TrialNbr))
            if size(table_normal_loading(idx_trial,:),1) ~= size(Ttemp,1)
                fprintf('Tables have different number of categorical values\n'); 
                fprintf('Advance carefully with new values\n'); 
            end
            table_normal_loading(idx_trial,:) = []; 
            table_normal_loading = [table_normal_loading;Ttemp];
        else
            table_normal_loading = [table_normal_loading;Ttemp];
        end
    end
end
fprintf(sprintf("data gather into table : [%s]\n",num2str(TRIAL_TARGET(trial_target_list))));
% to be saved now ! 
save(file_table_target,'table_normal_loading');
fprintf('data saved\n'); 


%% Visu segmentation 
trial_target = 1:Ntrial; 
list_force_target = [0.00,0.01,0.02,0.03,0.04,0.05,0.10,0.20,0.50,1:1:5];
% graph
fig_name = sprintf('visu segmentation %s',SUBJECT); 
fig = newfig(fig_name); Nrow = 4; Ncol = length(trial_target); 
fig.Position = [0    0.7938   50.8000   25.7440]; 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
ii = 0; 
for mm = trial_target
    ii = ii + 1; 
    
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>5,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    list_frame_num = 1:100; 
    frame_target = list_frame_num(1); 
    
    CAlast = T{mm}.ca{end}; 
    CAfit = myfitellipse(CAlast(:,1),CAlast(:,2)); 
    CAfit =[CAfit,repelem(mean(CAlast(:,3),'omitnan'),length(CAfit),1)];
    ca_idx_current = T{mm}.mask_in_border(:,frame_target);
    
    %ax(1)
    ax_idx = ii+0*Ncol;
    %centroid of ellipse at each frame
    m_ca = zeros(length(list_frame_num),3);
    for iframe = 1:length(list_frame_num)
        CA_current = T{mm}.ca{iframe}; 
        if size(CA_current,1) >= 10
            [CAfit_ii] = myfitellipse(CA_current(:,1),CA_current(:,2));
            plot(h_ax(ax_idx),CAfit_ii(:,1),CAfit_ii(:,2),'.-','MarkerSize',5);
            plot(h_ax(ax_idx),mean(CAfit_ii(:,1)),mean(CAfit_ii(:,2)),'k.-','MarkerSize',5);
            m_ca(iframe,:) = cat(2,mean(CAfit_ii),mean(CAlast(:,3)));
        end
    end
    text(h_ax(ax_idx),m_ca(end,1),min(CAfit_ii(:,2)),sprintf('trial #%d',T{mm}.trial_nbr),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',15);
    text(h_ax(ax_idx),min(CAfit_ii(:,1)),m_ca(end,2),sprintf('fitted ellipse\nat each frame'),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',15,'Rotation',90);
    validLogic_ca = m_ca(:,1) ~= 0; 
    plot(h_ax(ax_idx),m_ca(validLogic_ca,1),m_ca(validLogic_ca,2),'k-'); 
    scatter(h_ax(ax_idx),m_ca(validLogic_ca,1),m_ca(validLogic_ca,2),20,list_frame_num(validLogic_ca)','filled'); %list_frame_num(validLogic_ca)'
    p1 = [max(CAfit_ii(:,1)),max(CAfit_ii(:,2))];
    p2 = [m_ca(end,1),m_ca(end,2)];
    u = p2(1)-p1(1); v = p2(2)-p1(2);
    quiver(h_ax(ax_idx),p1(1),p1(2),u,v,'k','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','off');
    text(h_ax(ax_idx),p1(1),p1(2),sprintf('centroïds'),...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12,'Rotation',0);
    set(h_ax(ax_idx),'XLim',[-5,5],'YLim',[-5,15],...
        'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','View', [0 -90],'visible','off');
    colormap(h_ax(ax_idx),jet); 
    
    %ax(2)
    ax_idx = ii+1*Ncol; 
    h_rad = cell(1,Nrad_segmentation);
    PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
    scatter(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),1,[1 1 1].*0.5,"filled");
    for irad = 1:Nrad_segmentation
        h_rad{irad} = scatter(h_ax(ax_idx),PCcurrent(idx_IN_rad_all{mm}{irad},1),...
            PCcurrent(idx_IN_rad_all{mm}{irad},2),...
            5,color_set_array(irad,:),"filled");
    end
    plot(h_ax(ax_idx),m_ca(validLogic_ca,1),m_ca(validLogic_ca,2),'k-'); 
    scatter(h_ax(ax_idx),m_ca(validLogic_ca,1),m_ca(validLogic_ca,2),20,'k','filled'); %list_frame_num(validLogic_ca)'
%     legend(h_ax(ax_idx),[h_rad{:}],num2str((1:Nrad_segmentation)'));
    set(h_ax(ax_idx),'XLim',[-15,15],'YLim',[-5,15],...
        'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','View', [0 -90],'visible','off');
    set(h_ax(ax_idx),'clim',[1, list_frame_num(end)]); 
    text(h_ax(ax_idx),min(PCcurrent(:,1),[],'omitnan'),mean(PCcurrent(:,2),'omitnan'),sprintf('bottom view'),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',15,'Rotation',90);
    text(h_ax(ax_idx),min(PCcurrent(:,1),[],'omitnan')-12.5,max(PCcurrent(:,2),[],'omitnan'),sprintf('Fingertip segmentation'),...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Rotation',90);
    
    %ax(3)
    ax_idx = ii+2*Ncol;
    h_rad = cell(1,Nrad_segmentation);
    PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
    scatter3(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),PCcurrent(:,3),1,[1 1 1].*0.5,"filled");
    for irad = 1:Nrad_segmentation
        h_rad{irad} = scatter3(h_ax(ax_idx),PCcurrent(idx_IN_rad_all{mm}{irad},1),...
            PCcurrent(idx_IN_rad_all{mm}{irad},2),...
            PCcurrent(idx_IN_rad_all{mm}{irad},3),10,color_set_array(irad,:),"filled");
    end
    plot3(h_ax(ax_idx),m_ca(validLogic_ca,1),m_ca(validLogic_ca,2),m_ca(validLogic_ca,3),'k-'); 
    scatter3(h_ax(ax_idx),m_ca(validLogic_ca,1),m_ca(validLogic_ca,2),m_ca(validLogic_ca,3),20,list_frame_num(validLogic_ca)','filled'); 
%     legend(h_ax(ax_idx),[h_rad{:}],num2str((1:Nrad_segmentation)'));
    set(h_ax(ax_idx),'XLim',[-15,15],'YLim',[-5,15],...
        'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','View', [0 0],'visible','off');
    set(h_ax(ax_idx),'clim',[1, list_frame_num(end)]); 
    text(h_ax(ax_idx),min(PCcurrent(:,1),[],'omitnan'),mean(PCcurrent(:,3),'omitnan'),sprintf('front view'),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',15,'Rotation',90);

    
    %ax(4)
    ax_idx = ii+3*Ncol;
    logic_inside_rad_cur = logic_inside_rad_all{mm};
    bin_count_cum = sum(logic_inside_rad_cur,2);
%     bar(h_ax(ax_idx),list_frame_num,bin_count_cum(list_frame_num));
    pcolor(h_ax(ax_idx),~logic_inside_rad_cur(list_frame_num,:)')
    set(h_ax(ax_idx), 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
    colormap(h_ax(ax_idx),gray(2))
    grid(h_ax(ax_idx),'on');
    set(h_ax(ax_idx),'XLim',[list_frame_num(1),list_frame_num(end)],'YLim',[1 size(logic_inside_rad_cur,2)],'clim',[0 1]);
    h_ax(ax_idx).XLabel.String = 'frame#';
    h_ax(ax_idx).YLabel.String = 'Radial segmentation bin number';
    h_ax(ax_idx).Title.String = 'Inside contact Logic';
end
exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png')); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
face_target = 1:2; frame_target = 150;
ylim_set = [-10,10];
%figure settings
title_fig = sprintf("%s : test %s-%s-%s-%s-%dN-%1.2f-cum-fr:%d",SUBJECT,'speckles',...
    material,PHASE,NFCOND_SET,SPDDXLCOND_SET,frame_target); 
hf = newfig(title_fig);
set(hf,'Unit', 'Normalized','Position',[0    0.0278    1.0000    0.9009]);
set( hf, 'Color', 'White') ;
color_line_intersection = {'b','r','k','k','k','k','k'};  
% ax size
Nrow = 2; %Ncirc; 
Ncol = Ntrial; % boxplot total / map of regions
merge_element = {}; %mat2cell((1:2*length(face_target)),1,repelem(2,length(face_target))); 
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element,'shownum',0);

hold(h_ax,'on');
set(h_ax,'FontSize',12); 
itrial = 0; 
for mm = 1:Ntrial
    itrial = itrial +1; 
    ax_idx = mm; %+(xx-1)*Ncol;
    FCcurrent = FCmat_cum_last_all{mm}; 
    % z-line
    h_zline = yline(h_ax(ax_idx),0,'-','LineWidth',1);
    %boxplot 
    Ndata = size(FCcurrent,1); 
    legend_test = mat2cell(string(num2str((1:Nrad_segmentation)')),repelem(1,Nrad_segmentation),1); 
    x = mat2cell(FCcurrent,Ndata,repelem(Nrad_segmentation,Nface)); 
    x = x(face_target);
    h = boxplotGroup(h_ax(ax_idx),x,'groupLines',true,...
        'PrimaryLabels',FACEMEASURE(face_target), ...
        'SecondaryLabels',legend_test, ...
        'PlotStyle','Compact','GroupLabelType','Vertical',...
        'BoxStyle','filled',...
        'Colors',color_set_array,'GroupType','withinGroups');
    h_line_link = gobjects(1,length(face_target));
    for jj=1:length(face_target)
        medOuter = findobj(h.boxplotGroup(length(face_target)+1-jj),'tag','MedianOuter');
        xpos=flip([medOuter(:).XData]);
        h_line_link(jj)=plot(h_ax(ax_idx),xpos,median(x{jj},'omitnan'),'-',...
            'Color',color_line_intersection{jj},'LineWidth',2);
    end
    %Change to color of the median points. 
    set(h.xline,'LabelVerticalAlignment','bottom')
    medOuter = findobj(h.boxplotGroup,'tag','MedianOuter');
    set(medOuter,'MarkerFaceColor', 'k')
    medInner = findobj(h.boxplotGroup,'tag','MedianInner');
    set(medInner,'MarkerEdgeColor', 'w');
    %Change outlier symbol.
    outliers = findobj(h.boxplotGroup,'Tag','Outliers');
    set(outliers, 'Marker','.'); 
    %Change width of whiskers.
    set(findobj(h.boxplotGroup,'tag','Whisker'),'LineWidth',1)
    %Rotate xtick labels and font.
    h.axis.XTickLabelRotation = 90; 
%     h.axis.XAxis.FontName = 'fixedwidth';
    h.axis.XAxis.FontWeight = 'bold';
%     if xx~=Ncirc; set(h.axis.XAxis,'Visible','off'); end
    %Set the background color for each group.
    groupX = [cell2mat(get(h.xline,'value'));max(xlim(h.axis))];
    %Set YLim 
    h_ax(ax_idx).YLim = ylim_set; 
    yl = ylim(h.axis);
    ph = gobjects(1,numel(groupX)-1);
    for gg = 1:Nrad_segmentation-Nrad_out-1
        ph(gg) = patch(h.axis, ...
            groupX(gg+[0,1,1,0]),...
            yl([1,1,2,2]),...
            [1 1 1].*0.5,...
            'FaceAlpha',0.5,...
            'LineStyle','none');
    end
    %highlight the contour
    for gg = Nrad_segmentation-Nrad_out:Nrad_segmentation-Nrad_out
        ph(gg) = patch(h.axis, ...
            groupX(gg+[0,1,1,0]),...
            yl([1,1,2,2]),...
            'r',...
            'FaceAlpha',0.5,...
            'LineStyle','none');
    end
    for gg = Nrad_segmentation-Nrad_out+1:Nrad_segmentation
        ph(gg) = patch(h.axis, ...
            groupX(gg+[0,1,1,0]),...
            yl([1,1,2,2]),...
            [1 1 1],...
            'FaceAlpha',0.5,...
            'LineStyle','none');
    end
    
%     uistack(h_line_link,'bottom');
    uistack(h_zline,'bottom');
    uistack(ph,'bottom');
    %Turn on horizontal grid
    h.axis.YGrid = 'on'; 
    %Box off 
    box(h_ax(ax_idx),'off');
%     if (xx==1); title(h_ax(ax_idx),'Total (%)'); end
%     if xx==(floor(Ncirc/2)+1)
    h_ax(ax_idx).YLabel.String = 'Total (%)'; 
%     else
%         h_ax(ax_idx).YLabel.String = '';
%     end
    h_ax(ax_idx).YLabel.FontSize = 20; 
    h_ax(ax_idx).YLabel.FontWeight = 'bold'; 
    
%     if xx == 1 ; text(h_ax(ax_idx),mean(xpos)-0.5,ylim_set(2),sprintf('trial #%d',T{mm}.trial_nbr),...
%             'HorizontalAlignment','center','VerticalAlignment','bottom',...
%             'FontSize',12,'FontWeight','bold'); end
%     set(h_ax(ax_idx),'Color','None','visible','off');
%     set(findall(h_ax(ax_idx), 'type', 'text'),'visible','on')
    for iface = 1:length(face_target)
        ax_idx = itrial+1*Ncol;%iface*Ncol;  
        xdecalage = 20; 
        % map 
        PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
        PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
        FCcurrent = x{iface}; 
        ca_idx_current = T{mm}.mask_in_border(:,frame_target);
        
        for irad = 1:Nrad_segmentation
            ca_cur = findboundary(PCcurrent_segmented(:,:,irad),'connex_factor',0.0);
            plot(h_ax(ax_idx),ca_cur(:,1)+(iface-1)*xdecalage,ca_cur(:,2),'k-'); 
        end
        % graph 
        FCcurrent = sum(FCcurrent,2,'omitnan'); 
        FCcurrent(FCcurrent == 0) = NaN; 
        scatter(h_ax(ax_idx),PCcurrent(:,1)+(iface-1)*xdecalage,PCcurrent(:,2),5,[1 1 1].*0.7,"filled");
        scatter(h_ax(ax_idx),PCcurrent(:,1)+(iface-1)*xdecalage,PCcurrent(:,2),5,FCcurrent,"filled");
        colormap(h_ax(ax_idx),'jet');
        set(h_ax(ax_idx),'XLim',[-7,7+xdecalage],'YLim',[-5,20],'clim',[-10,10],...
            'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
    end
end
% allow out of plot data
set(h_ax,'Clipping','off');

%% Annulus strains 
fig_name = sprintf('%s - annulus rate and total strains',SUBJECT); 
fig = newfig(fig_name); 
set(fig,'Position',[3.7306   10.1865   13.4144   13.9435]); 
Nrow = 2; Ncol = 2; 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
itrial = 0; 
for mm = 1:Ntrial
    itrial = itrial + 1; 
    ylim_set = [-10,10];
    FCmat_annulus = FC_rate_val_annulus_all(:,:,mm);
    %ax(1)
    ax_idx = 1; 
    yline(h_ax(ax_idx),0,'k'); 
    plot(h_ax(ax_idx),t,squeeze(FCmat_annulus(:,1)),'Color',[1 1 1].*0.5,'LineWidth',2);
    plot(h_ax(ax_idx),t,squeeze(FCmat_annulus(:,2)),'Color',[1 1 1].*0.0,'LineWidth',2);
    set(h_ax(ax_idx),'XLim',[0,2],'YLim',ylim_set,'Clipping','off'); 
    set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off');
    ylabel(h_ax(ax_idx),'Rate (%/s)');
    mag_nf = NFCOND_SET*0.95;
    idx_nf = t(find(T{mm}.nf>mag_nf,1,'first')); 
    xline(h_ax(ax_idx),idx_nf);
    if mm == Ntrial
        text(h_ax(ax_idx),idx_nf,ylim_set(2),sprintf('%1.2fN',mag_nf),...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left'); 
    end
    
    %ax(2)
    ax_idx = 2; 
    yline(h_ax(ax_idx),0,'k'); 
    plot(h_ax(ax_idx),T{mm}.nf,squeeze(FCmat_annulus(:,1)),'Color',[1 1 1].*0.5,'LineWidth',2);
    plot(h_ax(ax_idx),T{mm}.nf,squeeze(FCmat_annulus(:,2)),'Color',[1 1 1].*0.0,'LineWidth',2);
    set(h_ax(ax_idx),'XLim',[0,NFCOND_SET],'YLim',ylim_set,'Clipping','off'); 
    set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off'); h_ax(ax_idx).YTickLabel = [];
    
    ylim_set = [-10,10];
    FCmat_annulus = FC_cum_val_annulus_all(:,:,mm);
    %ax(3)
    ax_idx = 3; 
    yline(h_ax(ax_idx),0,'k'); 
%     plot(h_ax(ax_idx),t,cumsum(FC_rate_val_annulus_all(:,1,mm)/50,'omitnan'));
    plot(h_ax(ax_idx),t,squeeze(FCmat_annulus(:,1)),'Color',[1 1 1].*0.5,'LineWidth',2);
    plot(h_ax(ax_idx),t,squeeze(FCmat_annulus(:,2)),'Color',[1 1 1].*0.0,'LineWidth',2);
    set(h_ax(ax_idx),'XLim',[0,2],'YLim',ylim_set,'Clipping','off'); 
    xlabel(h_ax(ax_idx),'Time (s)'); ylabel(h_ax(ax_idx),'Total (%)');
    xline(h_ax(ax_idx),idx_nf);
    if mm == Ntrial
        text(h_ax(ax_idx),idx_nf,ylim_set(2),sprintf('%1.2fN',mag_nf),...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left'); 
    end
    %ax(4)
    ax_idx = 4; 
    yline(h_ax(ax_idx),0,'k'); 
    plot(h_ax(ax_idx),T{mm}.nf,squeeze(FCmat_annulus(:,1)),'Color',[1 1 1].*0.5,'LineWidth',2);
    plot(h_ax(ax_idx),T{mm}.nf,squeeze(FCmat_annulus(:,2)),'Color',[1 1 1].*0.0,'LineWidth',2);
    set(h_ax(ax_idx),'XLim',[0,NFCOND_SET],'YLim',ylim_set,'Clipping','off');  h_ax(ax_idx).YTickLabel = [];
    xlabel(h_ax(ax_idx),'Force (N)'); 
    set(h_ax,'FontSize',12,'Color','None');
    
    
end
linkaxes([h_ax(1),h_ax(3)],'x'); 
linkaxes([h_ax(2),h_ax(4)],'x'); 


face_target = 1:2; 
trial_target = 1:Ntrial;
ptile_face = [1,99];
color_face = [[1 1 1].*0.5;...
              [1 1 1].*0.0];
          
facemaxrate = zeros(length(trial_target),length(list_frame_num),length(face_target)); 
idx_max_rate = zeros(length(trial_target),length(list_frame_num),length(face_target)); 
          
itrial = 0; 
for mm = trial_target
    itrial = itrial + 1; 
    list_frame_num = 1:T{mm}.Nframe;
    
    FCmat_cur = FCmat_rate{mm}; 
    for iface = 1:length(face_target)
        facemeasure_ii = FACEMEASURE{face_target(iface)}; 
        annulus_current = squeeze(FC_rate_val_annulus_all(:,face_target(iface),mm)); 
        frame_count = 0; 
        for iframe = list_frame_num 
            frame_count = frame_count + 1; 
            facemaxrate(itrial,frame_count,iface) = prctile(T{mm}.FC.(facemeasure_ii).rate(:,iframe),ptile_face(iface));
        end
    end
    
end

fig_name = sprintf('%s - annulus and max rate',SUBJECT); 
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.3380    0.3565    0.2641    0.4880]); 
Nrow = 1; Ncol = length(face_target); 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
ph = cell(1,length(face_target));
for iface = 1:length(face_target)
    ax_idx = iface;
    fmax_current = squeeze(facemaxrate(:,:,iface))';
    annulus_current = squeeze(FC_rate_val_annulus_all(:,face_target(iface),:));
    itrial = 0; 
    logic_valid = ~isnan(annulus_current(:)); 
    annulus_current(~logic_valid) = [];
    fmax_current(~logic_valid) = []; 
    xlim_set = [-40,40]; 
    xvec = linspace(xlim_set(1),xlim_set(2));
    
    %linear fit model
    lm = fitlm(fmax_current,annulus_current,'linear','RobustOpts','on');
    %plot
    ph{iface}=plot(h_ax(ax_idx),fmax_current,annulus_current,'.','Color','k');
    plot(h_ax(ax_idx),xvec,lm.Coefficients.Estimate(2)*xvec+lm.Coefficients.Estimate(1),'--','Color',[1 1 1].*0.5);
    
    set(h_ax(ax_idx),'XLim',xlim_set,'YLim',xlim_set,'clipping','off','DataAspectRatio',[1 1 1],'Color','None');
    ph_line = plot(h_ax(ax_idx),get(h_ax(ax_idx),'XLim'),get(h_ax(ax_idx),'YLim'),'-','Color','k');
    xlabel(h_ax(ax_idx),'Extremum');
    ylabel(h_ax(ax_idx),'Annulus');
    uistack(ph_line,'bottom');
    title(h_ax(ax_idx),FACEMEASURE{face_target(iface)}); 
end
% A tout temps, le maximum de epc1 est toujours très proche de la valeur
% moyenne de la bordure.

%%
trial_target = [1];
list_force_target = [0.05,0.1,0.2,0.5,1:5];
fig_name = sprintf('%s - annulus location with force',SUBJECT); 
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.0005    0.7417    0.9990    0.1796]); 
Nrow = trial_target; Ncol = length(list_force_target); 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
for mm = 1:length(trial_target)
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>5,1,'first')+1:end) = NaN;
    frame_target = find_closest_value(nf_format,list_force_target)';
    mask_in_annulus_current = T{mm}.mask_in_annulus(:,frame_target);
    mask_in_current = T{mm}.mask_in(:,frame_target);
    idx_frontal_target = T{mm}.mask_circ{1}(:,1)|T{mm}.mask_circ{2}(:,1)|T{mm}.mask_circ{3}(:,1)|T{mm}.mask_circ{4}(:,1)|T{mm}.mask_circ{5}(:,1);
    
    % map
    PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
    
    % graph
    xlim_map_set = [-10,10];
    ylim_map_set = [-10,15];
    
    for iframe = 1:length(list_force_target)
        %ax(iframe)
        PCcurrent_segmented = PCcurrent(mask_in_annulus_current(:,iframe)&idx_frontal_target,:,iframe);
        CA = T{mm}.ca{frame_target(iframe)};
        ax_idx = iframe+Ncol*(mm-1);
        scatter(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),1,[1 1 1].*0.8,"filled");
        scatter(h_ax(ax_idx),PCcurrent_segmented(:,1),PCcurrent_segmented(:,2),2,[0.8500 0.3250 0.0980],"filled");
        plot(h_ax(ax_idx),CA(:,1),CA(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
        set(h_ax(ax_idx),'XLim',xlim_map_set,'YLim',ylim_map_set,...
            'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
        text(h_ax(ax_idx),max(CA(:,1)),max(CA(:,2)),sprintf('at %1.2fN',list_force_target(iframe)),...
            'HorizontalAlignment','left','VerticalAlignment','middle');
    end
end


%% Face measure : multiple rad with force parameter 
% C = colormap(gray(round(Nrad_segmentation*(1.2)))); 
list_face_target_num = [1,2];
trial_target = 1:Ntrial; 
list_force_target = [0.5,1,2,4]; %[0.03,0.05,0.10,0.25,0.50,1:1:4,4.7];
ylim_set = [-10,10];

color_face = [[0 0 1];[0 1 0]];%[[1 1 1].*0.7;[1 1 1].*0.0]; %[1 1 1].*(linspace(0,0.7,Ntrial))';
fig_name = sprintf('%s - Epc1,Epc2 front - rate multiple rad trial 4',SUBJECT);%%1.3f,list_force_target); 
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.0734    0.2380    0.2474    0.6065]); 
Nrow = 7; Ncol = Nrad_segmentation;%length(list_force_target);%Ntrial; 
merge_element = {1:5*Ncol}; %{}; %
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');
itrial = 0; 
for mm = trial_target
    itrial = itrial + 1; 
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>NFCOND_SET,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    frame_target = list_frame_num(1); 
    
    FCmat_cur = FCmat_cum{mm};
    logic_inside_rad_cur = logic_inside_rad_all{mm};
%     FCmat_annulus = FC_cum_val_out_annulus_all(:,:,mm);
    ii = 0; 
%     yline(h_ax(ax_idx),0,'Color',[1 1 1].*0.5); 
    xvec = 1:Nrad_segmentation; 
    hp = cell(1,Nface); 
    ax_idx = 1;
    plot(h_ax(ax_idx),repelem(1:Nrad_segmentation,2,1),[ylim_set(1).*ones(1,Nrad_segmentation);...
        satur(min(FCmat_cur(list_frame_num,1:Nrad_segmentation),[],1),'level',ylim_set(1),'method','low')],'-','Color',[1 1 1].*0.5);
    for iframe = list_frame_num
        ii = ii + 1;  
        %itrial; 
        hpline = yline(h_ax(ax_idx),0);
        logic_inside_current = logic_inside_rad_cur(iframe,:);
        FCmat_cur_in = FCmat_cur(iframe,:); FCmat_cur_out = FCmat_cur(iframe,:);
        FCmat_cur_in(repmat(~logic_inside_current,1,Nface)) = NaN; 
        FCmat_cur_out(repmat(logic_inside_current,1,Nface)) = NaN; 
        
        kk = 0; 
        for iface = list_face_target_num
            kk = kk + 1; 
            %face measure ii
            hp{kk} = plot(h_ax(ax_idx),xvec,FCmat_cur(iframe,(1:Nrad_segmentation)+(iface-1)*Nrad_segmentation)','-','LineWidth',2,'Color',[color_face(iface,:),0.2]); %[1 1 1].*color_scale_grey(ii));
%             scatter(h_ax(ax_idx),xvec,FCmat_cur_in((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation),30,repelem(nf_format(iframe),Nrad_segmentation),'filled');
%             scatter(h_ax(ax_idx),xvec,FCmat_cur_out((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation),30,repelem(nf_format(iframe),Nrad_segmentation),'MarkerFaceColor','w','LineWidth',1);
            scatter(h_ax(ax_idx),xvec,FCmat_cur_in((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation),30,'k','filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
            scatter(h_ax(ax_idx),xvec,FCmat_cur_out((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation),30,'k','MarkerFaceColor','w','LineWidth',1,'MarkerEdgeAlpha',0.2);
           	if mm == trial_target(end)
            text(h_ax(ax_idx),Nrad_segmentation*1.01,FCmat_cur(iframe,iface*Nrad_segmentation),...
                sprintf('%1.1fN',list_force_target(ii)),...
                'HorizontalAlignment','left',...
                'FontSize',8);
            end
        end
%         if itrial == 1
%         text(h_ax(ax_idx),round(Nrad_segmentation/2),35,sprintf('%1.1fN',list_force_target(ii)),...
%             'FontSize',12,...
%             'HorizontalAlignment','center');
%         end
        uistack([hpline],'bottom');
        set(h_ax(ax_idx),'XLim',[0,Nrad_segmentation+1],'YLim',ylim_set,'Clipping','off');%[-0.025,0.025][-1.5,1.5][0,2e-3]
        ylim_set = get(h_ax(ax_idx),'YLim');
        %     for kk = 1:Nrad_segmentation
        %         plot(h_ax(ax_idx),[kk,kk],[ylim_set(1),FCmat_cur(end,kk)],'k-','LineWidth',1);
        %     end
        colormap(h_ax(ax_idx),'jet'); %[1 1 1].*color_scale_grey');
        if false%ax_idx == length(h_ax)
            h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');
            h_cb.Title.String=sprintf('Force (N)');
            
        end
        if ax_idx ~= 1
            set(get(h_ax(ax_idx), 'YAxis'), 'Visible', 'off');
        end
        set(h_ax(ax_idx),'YGrid','on'); 
        set(h_ax(ax_idx),'clim',[0,NFCOND_SET],'Color','none');
        h_ax(ax_idx).XLabel.String = 'Meridional Bin';
        h_ax(ax_idx).YLabel.String = 'Strain (%)';
        %     h_ax(ax_idx).XTickLabel = [];
        set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off');
%         lgd_label = cat(2,FACEMEASURE(list_face_target_num),{'in','out'});
%         qw = cell(1,2);
%         qw{1} = plot([NaN],'ko','MarkerSize',10,'MarkerFaceColor','k');
%         qw{2} = plot([NaN],'ko','MarkerSize',10,'MarkerFaceColor','w');
%         legend(h_ax(ax_idx),[hp{:},qw{:}],lgd_label,'Location','northwest','box','off','NumColumns',2);
    
    end
    if itrial == 1
        % map 
        frame_target = 150; 
        PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
        PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
        ca_idx_current = T{mm}.mask_in_border(:,frame_target);

        % graph 
        for irad = 1:Nrad_segmentation
            ax_idx = 5*Ncol+irad;
            if mod(ax_idx,2)==0 
                ax_idx = ax_idx+Ncol;
            end
            scatter3(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),PCcurrent(:,3),1,[1 1 1].*0.8,"filled");
            scatter3(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),PCcurrent_segmented(:,3,irad),5,[0.8500 0.3250 0.0980],"filled");
    %         plot(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
            set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
                'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
            colormap(h_ax(ax_idx),'jet');
        end

        set(h_ax,'visible','off');
        set(h_ax(1),'visible','on');
        set(h_ax,'FontSize',12); 
    end
end
set(h_ax,'FontSize',12); 
exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png')); 

%%
list_face_target_num = [1,2];
list_force_target = [0.05,0.1,0.2,0.5,1,2,4];
nanvec = NaN(Nrad_segmentation,length(list_force_target),length(list_face_target_num),Ntrial); 
FC_all = nanvec; FC_all_in = nanvec; FC_all_out = nanvec;
itrial = 0; 
for mm = 1%:Ntrial 
    itrial = itrial + 1; 
    FCmat_cur = FCmat_cum{mm};
    logic_inside_rad_cur = logic_inside_rad_all{mm};
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>NFCOND_SET,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)'
    for iframe = 1:length(list_frame_num)
        logic_inside_current = logic_inside_rad_cur(list_frame_num(iframe),:);
        FCmat_cur_in = FCmat_cur(list_frame_num(iframe),:); 
        FCmat_cur_out = FCmat_cur(list_frame_num(iframe),:);
        FCmat_cur_in(repmat(~logic_inside_current,1,Nface)) = NaN; 
        FCmat_cur_out(repmat(logic_inside_current,1,Nface)) = NaN; 
        for iface = list_face_target_num
            kk = kk + 1; 
            FC_all(:,iframe,iface,itrial) = FCmat_cur(list_frame_num(iframe),(1:Nrad_segmentation)+(iface-1)*Nrad_segmentation)';
            FC_all_in(:,iframe,iface,itrial) = FCmat_cur_in((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation)';
            FC_all_out(:,iframe,iface,itrial) = FCmat_cur_out((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation)';
        end
    end
end

ylim_set = [-10,10];

color_face = [[0 0 1];[0 1 0]];%[[1 1 1].*0.7;[1 1 1].*0.0]; %[1 1 1].*(linspace(0,0.7,Ntrial))';
fig_name = sprintf('%s - Epc1,Epc2 front - rate multiple rad mean %1.2f',SUBJECT,list_force_target(end));%%1.3f,list_force_target); 
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.0734    0.2380    0.2474    0.6065]); 
Nrow = 7; Ncol = Nrad_segmentation;%length(list_force_target);%Ntrial; 
merge_element = {1:5*Ncol}; %{}; %
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');

xvec = 1:Nrad_segmentation; 
hp = cell(1,Nface); 

ax_idx = 1;
hpline = yline(h_ax(ax_idx),0);
kk = 0; 
for iface = list_face_target_num
    kk = kk + 1; 
    %face measure
    if length(list_frame_num) ~=1
        FC_current = mean(squeeze(FC_all(:,:,kk,:)),3);
        FC_all_in_current = mean(squeeze(FC_all_in(:,:,kk,:)),3);
        FC_all_out_current = mean(squeeze(FC_all_out(:,:,kk,:)),3);
    else
        FC_current = mean(squeeze(FC_all(:,:,kk,:)),2);
        FC_all_in_current = mean(squeeze(FC_all_in(:,:,kk,:)),2);
        FC_all_out_current = mean(squeeze(FC_all_out(:,:,kk,:)),2);
    end
    if iface == 1
        plot(h_ax(ax_idx),repelem(1:Nrad_segmentation,2,1),[ylim_set(1).*ones(1,Nrad_segmentation);...
            satur((min(FC_current(:,1:length(list_frame_num)),[],2))','level',ylim_set(1),'method','low')],'-','Color',[1 1 1].*0.5);
    end
    hp{kk} = plot(h_ax(ax_idx),xvec,FC_current,'-','LineWidth',2,'Color',[color_face(kk,:)]); 
    scatter(h_ax(ax_idx),xvec,FC_all_in_current,30,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
    scatter(h_ax(ax_idx),xvec,FC_all_out_current,30,'k','MarkerFaceColor','w','LineWidth',1,'MarkerEdgeAlpha',1);
    for iframe = 1:length(list_frame_num)
%         text(h_ax(ax_idx),Nrad_segmentation*1.01,FC_current(end,iframe),...
%             sprintf('%1.2fN',list_force_target(iframe)),...
%             'HorizontalAlignment','left',...
%             'FontSize',8);
    end
end
uistack([hpline],'bottom');
set(h_ax(ax_idx),'XLim',[0,Nrad_segmentation+1],'YLim',ylim_set,'Clipping','off');%[-0.025,0.025][-1.5,1.5][0,2e-3]
ylim_set = get(h_ax(ax_idx),'YLim');

colormap(h_ax(ax_idx),'jet'); 
if ax_idx ~= 1
    set(get(h_ax(ax_idx), 'YAxis'), 'Visible', 'off');
end
set(h_ax(ax_idx),'YGrid','on'); 
set(h_ax(ax_idx),'clim',[0,NFCOND_SET],'Color','none');
h_ax(ax_idx).XLabel.String = 'Meridional Bin';
h_ax(ax_idx).YLabel.String = 'Strain (%)';
set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off');

% map
mm = 1; 
frame_target = 150;
PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
ca_idx_current = T{mm}.mask_in_border(:,frame_target);

% graph
for irad = 1:Nrad_segmentation
    ax_idx = 5*Ncol+irad;
    if mod(ax_idx,2)==0
        ax_idx = ax_idx+Ncol;
    end
    scatter3(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),PCcurrent(:,3),1,[1 1 1].*0.8,"filled");
    scatter3(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),PCcurrent_segmented(:,3,irad),5,[0.8500 0.3250 0.0980],"filled");
    set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
        'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
    colormap(h_ax(ax_idx),'jet');
end

set(h_ax,'visible','off');
set(h_ax(1),'visible','on');
set(h_ax,'FontSize',12); 
exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png')); 

%%
ylim_set = [-10,10];
list_force_target = [0.5,1,5];
color_face = [[1 1 1].*0.7;[1 1 1].*0.0]; %[1 1 1].*(linspace(0,0.7,Ntrial))';
fig_name = sprintf('%s - Epc1,Epc2 front - cum multiple rad last',SUBJECT);%%1.3f,list_force_target); 
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.3380    0.3565    0.2641    0.4880]); 
Nrow = 7; Ncol = Nrad_segmentation;%Ntrial; 
merge_element = {1:5*Ncol}; %{}; %
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');
itrial = 0; 
for mm = 1:Ntrial
    itrial = itrial + 1; 
    
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>NFCOND_SET,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    frame_target = list_frame_num(1); 
    
    ax_idx = 1;%itrial; 
    FCmat_cur = FCmat_cum{mm};
    logic_inside_rad_cur = logic_inside_rad_all{mm};
%     FCmat_annulus = FC_cum_val_out_annulus_all(:,:,mm);
    ii = 0; 
    yline(h_ax(ax_idx),0,'Color',[1 1 1].*0.5); 
    xvec = 1:Nrad_segmentation; 
    hp = cell(1,Nface); 
    plot(h_ax(ax_idx),repelem(1:Nrad_segmentation,2,1),[ylim_set(1).*ones(1,Nrad_segmentation);...
        satur(min(FCmat_cur(list_frame_num,1:Nrad_segmentation),[],1),'level',ylim_set(1),'method','low')],'-','Color',[1 1 1].*0.5);
    for iframe = list_frame_num
        ii = ii + 1;  
        logic_inside_current = logic_inside_rad_cur(iframe,:);
        FCmat_cur_in = FCmat_cur(iframe,:); FCmat_cur_out = FCmat_cur(iframe,:);
        FCmat_cur_in(repmat(~logic_inside_current,1,Nface)) = NaN; 
        FCmat_cur_out(repmat(logic_inside_current,1,Nface)) = NaN; 
        kk = 0; 
        for iface = list_face_target_num
            kk = kk + 1; 
            %face measure ii
%             hp{kk} = plot(h_ax(ax_idx),xvec,FCmat_cur(iframe,(1:Nrad_segmentation)+(iface-1)*Nrad_segmentation)','-','LineWidth',2,'Color',color_face(iface,:)); %[1 1 1].*color_scale_grey(ii));
%             scatter(h_ax(ax_idx),xvec,FCmat_cur_in((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation),30,repelem(nf_format(iframe),Nrad_segmentation),'filled');
%             scatter(h_ax(ax_idx),xvec,FCmat_cur_out((1:Nrad_segmentation)+(iface-1)*Nrad_segmentation),30,repelem(nf_format(iframe),Nrad_segmentation),'MarkerFaceColor','w','LineWidth',1);
%            	text(h_ax(ax_idx),Nrad_segmentation*1.01,FCmat_cur(iframe,iface*Nrad_segmentation),...
%                 sprintf('%1.2f',nf_format(iframe)),...
%                 'HorizontalAlignment','left',...
%                 'FontSize',8);
        end
    end
    set(h_ax(ax_idx),'XLim',[0,Nrad_segmentation+1],'YLim',ylim_set,'Clipping','off');%[-0.025,0.025][-1.5,1.5][0,2e-3]
    ylim_set = get(h_ax(ax_idx),'YLim');
%     for kk = 1:Nrad_segmentation
%         plot(h_ax(ax_idx),[kk,kk],[ylim_set(1),FCmat_cur(end,kk)],'k-','LineWidth',1); 
%     end
    colormap(h_ax(ax_idx),'jet'); %[1 1 1].*color_scale_grey'); 
%     h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');
%     h_cb.Title.String=sprintf('Force (N)'); 
    set(h_ax(ax_idx),'clim',[0,NFCOND_SET],'Color','none'); 
    h_ax(ax_idx).XLabel.String = 'Meridional Bin';  
    h_ax(ax_idx).YLabel.String = 'Rate (%/s)';
%     h_ax(ax_idx).XTickLabel = [];
    set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off');
    if itrial == 1
        % map 
        frame_target = 150; 
        PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
        PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
        ca_idx_current = T{mm}.mask_in_border(:,frame_target);

        % graph 
        for irad = 1:Nrad_segmentation
            ax_idx = 5*Ncol+irad;
            if mod(ax_idx,2)==0 
                ax_idx = ax_idx+Ncol;
            end
            scatter3(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),PCcurrent(:,3),1,[1 1 1].*0.8,"filled");
            scatter3(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),PCcurrent_segmented(:,3,irad),2,[0.8500 0.3250 0.0980],"filled");
    %         plot(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
            set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
                'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
            colormap(h_ax(ax_idx),'jet');
        end

        set(h_ax,'visible','off');
        set(h_ax(1),'visible','on');
        set(h_ax,'FontSize',12); 
    end
end
%mean val 
for iframe = list_frame_num
    kk = 0; 
for iface = list_face_target_num
kk = kk+1; 
val = zeros(Ntrial,Nrad_segmentation);
logic_inside_mean_all = NaN(Ntrial,Nrad_segmentation); 
for mm = 1:Ntrial
    val(mm,:) = FCmat_cum{mm}(iframe,(1:Nrad_segmentation)+(iface-1)*Nrad_segmentation);
    plot(h_ax(1),xvec,val(mm,:)','-','LineWidth',1,'Color',[1 1 1].*0.7); %[1 1 1].*color_scale_grey(ii));
    logic_inside_mean_all(mm,:) = logic_inside_rad_all{mm}(iframe,:);
end 
meanval = mean(val); 
ml = find(gradient(sum(logic_inside_mean_all,1)>2)~=0,1,'first'); 
hp{kk} = plot(h_ax(1),xvec,meanval','-','LineWidth',3,'Color',color_face(iface,:)); %[1 1 1].*color_scale_grey(ii));

plot(h_ax(1),xvec(ml),meanval(ml),'.','MarkerSize',20,'Color',[0.8500 0.3250 0.0980])
end
end

lgd_label = cat(2,FACEMEASURE(list_face_target_num));
legend(h_ax(ax_idx),[hp{:}],lgd_label,'Location','northwest','box','off','NumColumns',1); 
    

exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png')); 


%% Face measure : multiple force with rad parameter 
% C = colormap(gray(round(Nrad_segmentation*(1.2)))); 
frame_target = 1; 
list_face_target_num = [1,2];
list_force_target = [0,0.5,1:5];
list_rad_num = 1:Nrad_in; 
color_face = [[1 1 1].*0.7;...
              [1 1 1].*0.0];
ylim_set = [-10,10]; 

fig_name = sprintf('%s - Curv north - cum multiple force with rad',SUBJECT); 
fig = newfig(fig_name); 
set(fig,'Position',[3.7306   10.1865   13.4144   13.9435]); 
Nrow = 8; Ncol = Nrad_segmentation; 
merge_element = {1:5*Ncol};
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');
itrial = 0; 
for mm = 1%:Ntrial
    itrial = itrial + 1; 
    ax_idx = 1; 
    FCmat_cur = FCmat_cum{mm};
    logic_inside_rad_cur = logic_inside_rad_all{mm};
%     FCmat_annulus = FC_cum_val_out_annulus_all(:,:,mm);
    yline(h_ax(ax_idx),0,'Color',[1 1 1].*0.5); 
    
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>NFCOND_SET,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    
    xvec_nf = nf_all(list_frame_num,mm); 
    xvec = xvec_nf; 
    hp = cell(1,Nface); 
    
    ii = 0; 
    for irad = list_rad_num
        ii = ii + 1;  
        kk = 0; 
        for iface = list_face_target_num
            kk = kk + 1; 
            logic_inside_current = logic_inside_rad_cur(list_frame_num,irad);
            FCmat_current = FCmat_cur(list_frame_num,irad+(iface-1)*Nrad_segmentation);
            FCmat_cur_in = FCmat_current; FCmat_cur_out = FCmat_current;
            FCmat_cur_in(~logic_inside_current) = NaN; 
            FCmat_cur_out(logic_inside_current) = NaN; 
            %face measure ii
            hp{kk} = plot(h_ax(ax_idx),xvec,FCmat_current,'-','LineWidth',2,'Color',color_face(kk,:)); %[1 1 1].*color_scale_grey(ii));
            scatter(h_ax(ax_idx),xvec,FCmat_cur_in,30,repelem(irad,1*length(xvec),1),'filled');
            scatter(h_ax(ax_idx),xvec,FCmat_cur_out,30,repelem(irad,1*length(xvec),1),'Marker','s','MarkerFaceColor','w','LineWidth',2);
           	text(h_ax(ax_idx),xvec(end)+0.1,FCmat_current(end),...
                sprintf('%d',list_rad_num(ii)),...
                'HorizontalAlignment','left',...
                'FontSize',8);
        end
    end
    set(h_ax(ax_idx),'XLim',[0,xvec_nf(end)+0.5],'YLim',ylim_set,'Clipping','off');%[-0.025,0.025][-1.5,1.5]
    ylim_set = get(h_ax(ax_idx),'YLim');
%     for kk = 1:Nrad_segmentation
%         plot(h_ax(ax_idx),[kk,kk],[ylim_set(1),FCmat_cur(end,kk)],'k-','LineWidth',1); 
%     end
    cmap = flip(colormap(summer(length(list_rad_num))));
    colormap(h_ax(ax_idx),cmap); %[1 1 1].*color_scale_grey'); 
    h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');
    h_cb.Title.String=sprintf('Meridional\n Bin (-)'); 
    set(h_ax(ax_idx),'clim',[0,length(list_rad_num)],'Color','none'); 
    h_ax(ax_idx).XLabel.String = 'Force (N)';  
    h_ax(ax_idx).YLabel.String = "Total (%)";%,FACEMEASURE{list_face_target_num(1)});
%     h_ax(ax_idx).XTickLabel = [];
%     set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off');
    lgd_label = cat(2,FACEMEASURE(list_face_target_num),{'in','out'});
    qw{1} = plot([NaN],'ko','MarkerSize',10,'MarkerFaceColor','k'); 
    qw{2} = plot([NaN],'ko','MarkerSize',10,'MarkerFaceColor','w'); 
    legend(h_ax(ax_idx),[hp{:},qw{:}],lgd_label,'Location','northwest','box','off','NumColumns',2); 
  
    if itrial == 1
        % map 
        PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
        PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
        ca_idx_current = T{mm}.mask_in_border(:,frame_target);

        % graph 
        for irad = list_rad_num
            ax_idx = 6*Ncol+irad;
            if mod(ax_idx,2)==0 
                ax_idx = ax_idx+Ncol;
            end
            scatter(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),1,[1 1 1].*0.8,"filled");
            scatter(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),2,[0.8500 0.3250 0.0980],"filled");
    %         plot(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
            set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
                'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
            colormap(h_ax(ax_idx),'jet');
        end

        set(h_ax,'visible','off');
        set(h_ax(1),'visible','on');
        set(h_ax,'FontSize',12); 
    end
end
exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png'),...
                'BackgroundColor','none'); 


%% Profiles with forces

list_force_target = [0.01,0.015,0.02,0.03,0.05,0.10,0.20,0.50,1:1:5];
iface_target = 1; 
fig_name = sprintf('%s - profile multiple rad with force %1.3f',SUBJECT,list_force_target); 
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.0734    0.3565    0.2641    0.4880]); 
Nrow = 7; Ncol = Nrad_segmentation; 
merge_element = {1:5*Ncol};
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');
itrial = 0; 
for mm = 1%:Ntrial
    itrial = itrial + 1; 
    
    
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>5,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    frame_target = list_frame_num(1); 
    
    ax_idx = 1; 
    FCmat_cur = squeeze(FCmat_xyz{mm}(:,3,:))+1.5;
    logic_inside_rad_cur = logic_inside_rad_all{mm};
%     FCmat_annulus = FC_cum_val_out_annulus_all(:,:,mm);
    color_scale_grey = 1-linspace(0.10,1.00,(73-12)/10+1);
    ii = 0; 
    yline(h_ax(ax_idx),0,'Color',[1 1 1].*0.5); 
    x_vec = 1:Nrad_segmentation; 
    for iframe = list_frame_num
        ii = ii + 1;  
        logic_inside_current = logic_inside_rad_cur(iframe,:);
        FCmat_cur_in = FCmat_cur(iframe,:); FCmat_cur_out = FCmat_cur(iframe,:);
        FCmat_cur_in(~logic_inside_current) = NaN; 
        FCmat_cur_out(logic_inside_current) = NaN; 
        hp1 = plot(h_ax(ax_idx),1:Nrad_segmentation,FCmat_cur(iframe,1:Nrad_segmentation)','-','LineWidth',1,'Color',[1 1 1].*0.7); %[1 1 1].*color_scale_grey(ii));
        scatter(h_ax(ax_idx),x_vec,FCmat_cur_in,30,FCmat_rate{mm}(iframe,x_vec+(iface_target-1)*Nrad_segmentation),'filled');
        scatter(h_ax(ax_idx),x_vec,FCmat_cur_out,30,FCmat_rate{mm}(iframe,x_vec+(iface_target-1)*Nrad_segmentation),'LineWidth',2,'MarkerFaceColor','w');
        text(h_ax(ax_idx),Nrad_segmentation,FCmat_cur(iframe,Nrad_segmentation),...
            sprintf('%1.3f',list_force_target(ii)),...
            'HorizontalAlignment','left',...
            'FontSize',8);
    end
    set(h_ax(ax_idx),'XLim',[0,Nrad_segmentation+1],'YLim',[0,5],'Clipping','off');
    ylim_set = get(h_ax(ax_idx),'YLim');
    for kk = 1:Nrad_segmentation
        plot(h_ax(ax_idx),[kk,kk],[ylim_set(1),FCmat_cur(end,kk)],'k-','LineWidth',1); 
    end
    colormap(h_ax(ax_idx),'jet'); %[1 1 1].*color_scale_grey'); 
    h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');
    h_cb.Title.String=sprintf('Curvature (-)'); 
    set(h_ax(ax_idx),'clim',[-10,10],'Color','none'); 
    h_ax(ax_idx).XLabel.String = 'Meridional Bin';  
    h_ax(ax_idx).YLabel.String = 'Z-axis (mm)';
%     h_ax(ax_idx).XTickLabel = [];
    set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off');
%     legend(h_ax(ax_idx),[hp1],{'Profile'},'Location','northwest','box','off'); 
    
    
    
    % map 
    if itrial == 1
        PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
        PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
        ca_idx_current = T{mm}.mask_in_border(:,frame_target);

        % graph 
        for irad = 1:Nrad_segmentation
            ax_idx = 5*Ncol+irad;
            if mod(ax_idx,2)==0 
                ax_idx = ax_idx+Ncol;
            end
            scatter(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),1,[1 1 1].*0.8,"filled");
            scatter(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),2,[0.8500 0.3250 0.0980],"filled");
    %         plot(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
            set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
                'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
            colormap(h_ax(ax_idx),'jet');
        end

        set(h_ax,'visible','off');
        set(h_ax(1),'visible','on');
        set(h_ax,'FontSize',12); 
    end
end
exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png'),...
                'BackgroundColor','none'); 
%% EpcX vs curvature 
list_face_target_num = 1; 
list_frame_num = 1:Nframe; %27:5:97;
list_rad_num = 1:Nrad_in; %+1:Nrad_segmentation;%+1:Nrad_in+Nrad_out];
fig_name = sprintf('Epc1 vs curvature %s',SUBJECT); 
fig = newfig(fig_name); 

set(fig,'Unit','Normalized','Position',[0.0734    0.4133    0.2641    0.4392]); 
Nrow = 1; Ncol = length(list_face_target_num); 
h_ax = subplot_ax(Nrow,Ncol); hold('on'); 
for mm = 1%:Ntrial
    itrial = itrial + 1; 
    FCmat_cur = FCmat_cum{mm};
    ii = 0; 
    for irad = list_rad_num
        ii = ii + 1;  
        hp = cell(1,length(list_face_target_num));
        kk = 0; 
        face_count = 0; 
        for iface = list_face_target_num
            face_count = face_count + 1; 
            ax_idx = face_count;  
            kk = kk + 1; 
            logic_inside_current = logic_inside_rad_cur(list_frame_num,irad);
            FCmat_curv_current = FCmat_cum{mm}(list_frame_num,irad+(4-1)*Nrad_segmentation);
            FCmat_curv_cur_in = FCmat_curv_current; FCmat_curv_cur_out = FCmat_curv_current;
            FCmat_curv_cur_in(~logic_inside_current) = NaN; 
            FCmat_curv_cur_out(logic_inside_current) = NaN; 
            xvec = FCmat_curv_current;
            %face measure ii
            hp{kk} = plot(h_ax(ax_idx),T{mm}.nf,FCmat_curv_current,'-','LineWidth',1,'Color',color_face(kk,:)); 
            scatter(h_ax(ax_idx),T{mm}.nf,FCmat_curv_cur_in,30,FCmat_rate{mm}(list_frame_num,irad+(iface-1)*Nrad_segmentation),'filled');
            scatter(h_ax(ax_idx),T{mm}.nf,FCmat_curv_cur_out,30,FCmat_rate{mm}(list_frame_num,irad+(iface-1)*Nrad_segmentation),'Marker','s','MarkerFaceColor','w','LineWidth',2);
%             scatter(h_ax(ax_idx),T{mm}.nf,FCmat_curv_current,30,FCmat_rate{mm}(list_frame_num,irad+(iface-1)*Nrad_segmentation),'filled');

            h_ax(ax_idx).XLim = [0,NFCOND_SET]; 
            h_ax(ax_idx).YLim = [0,3e-3];
            set(h_ax(ax_idx),'clim',[-15,15],'Color','none','clipping','off'); 
            h_ax(ax_idx).XLabel.String = 'Force (N)'; 
            h_ax(ax_idx).YLabel.String = 'Curvature (-)'; 
            h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');
            h_cb.Title.String = sprintf('Rate (%/s)'); 
            cmap = colormap(jet(30));
            colormap(h_ax(ax_idx),cmap); 
            h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');        
        end
    end
end


%% Contact Area Evolution 
list_force_target = [0.05,0.1,0.2,0.5,1,2,4];%[0.02,0.03,0.04,0.05,0.10,0.20,0.50,1:1:5];
A = zeros(1,Nframe); 

fig_name = sprintf('Contact Area Evolution %s',SUBJECT);%,list_force_target);  %1.3fN
fig = newfig(fig_name); 
set(fig,'Unit','Normalized','Position',[0.1208    0.0287    0.1672    0.8926/3]); 
Nrow = 1; Ncol = 1; 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
for mm = 1:Ntrial
    for iframe = 1:Nframe
        CA_current = T{mm}.ca{iframe};
        if length(CA_current) > 3
            CAfit = myfitellipse(CA_current(:,1),CA_current(:,2)); 
            [x,y] = poly2cw(CAfit(:,1),CAfit(:,2));
            A(iframe) = polyarea(x,y);
            if abs(A(iframe))>500
                A(iframe) = NaN; 
            end
        end
    end
    
    nf_format = T{mm}.nf;
    nf_format = nf_format(1:find(nf_format>NFCOND_SET,1,'first'));
    list_frame_num = find_closest_value(nf_format,list_force_target);
    
%     displ = squeeze(FC_cum_val_all(:,3,1,4,mm));
%     displ = displ;% - displ(idx_touch_all(mm)); 
    
%     ax_idx = 1; 
%     plot(h_ax(ax_idx),t,displ,'k.-','MarkerSize',8);
%     ylim(h_ax(ax_idx),[2.45,2.7]); 
%     set(h_ax(ax_idx),'clipping','on');
%     h_ax(ax_idx).XLabel.String = 'Time [s]';
%     h_ax(ax_idx).YLabel.String = 'Indentation [mm]';
%     ax_idx = 1; 
%     plot(h_ax(ax_idx),t(1:length(nf_format)),nf_format,'k.-','MarkerSize',8);
%     plot(h_ax(ax_idx),t(list_frame_num),nf_all(list_frame_num,mm),'ro','MarkerSize',5,'MarkerFaceColor','r');
%     ylim(h_ax(ax_idx),[0,NFCOND_SET]); set(h_ax(ax_idx),'clipping','off');
%     h_ax(ax_idx).XLabel.String = 'Time [s]';
%     h_ax(ax_idx).YLabel.String = 'Force [N]';
%     ax_idx = 2; 
%     plot(h_ax(ax_idx),t,A,'k.-','MarkerSize',8);
%     plot(h_ax(ax_idx),t(list_frame_num),A(list_frame_num),'ro','MarkerSize',5,'MarkerFaceColor','r');
%     h_ax(ax_idx).XLabel.String = 'Time [s]';
%     h_ax(ax_idx).YLabel.String = 'Area [mm2]';
%     xlim(h_ax(ax_idx),[0,2]); 
    ax_idx = 1; 
    plot(h_ax(ax_idx),nf_format,A(1:length(nf_format)),'k.-','MarkerSize',8);
    plot(h_ax(ax_idx),nf_all(list_frame_num,mm),A(list_frame_num),'ro','MarkerSize',5,'MarkerFaceColor','r');
    xlim(h_ax(ax_idx),[0,NFCOND_SET]); set(h_ax(ax_idx),'clipping','off'); 
    h_ax(ax_idx).XLabel.String = 'Force [N]';
    h_ax(ax_idx).YLabel.String = 'Area [mm2]';
%     ax_idx = 6; 
%     plot(h_ax(ax_idx),displ,A,'k.-','MarkerSize',8);
%     xlim(h_ax(ax_idx),[2.45,2.7]); set(h_ax(ax_idx),'clipping','on'); 
%     h_ax(ax_idx).XLabel.String = 'Indentation [mm]';
%     h_ax(ax_idx).YLabel.String = 'Area [mm2]';
%     ax_idx = 2; 
%     plot(h_ax(ax_idx),displ,T{mm}.nf,'k.-','MarkerSize',8);
%     xlim(h_ax(ax_idx),[2.45,2.70]); 
%     set(h_ax(ax_idx),'clipping','on'); 
%     h_ax(ax_idx).XLabel.String = 'Indentation [mm]';
%     h_ax(ax_idx).YLabel.String = 'Force [N]';
    
end
grid(h_ax,'on');
set(h_ax,'FontSize',12,'Color','none');


exportgraphics(fig,sprintf('fig/%s.%s',fig_name,'png')); 




%%
frame_target = 86; 
for mm = 1%:Ntrial 
    PClast = T{mm}.FC.xyz(:,:,end);
    CAlast = T{mm}.ca{end}; 
    idx_frontal_last = T{mm}.mask_circ{1}(:,1)|T{mm}.mask_circ{2}(:,1)|T{mm}.mask_circ{3}(:,1)|T{mm}.mask_circ{4}(:,1)|T{mm}.mask_circ{5}(:,1); 
    % fit ellipse from last contact contour
    [CAfit,centerEl] = myfitellipse(CAlast(:,1),CAlast(:,2)); 
    %
    CAfit_rad = cell(1,Nrad_segmentation);
    idx_IN_rad = cell(1,Nrad_segmentation);
    idx_IN_allprevious = false(length(PClast),1);
    for irad = 1:Nrad_segmentation
        CAfit_rad{irad} = (CAfit - centerEl)*(irad/Nrad_in)+centerEl; 
        idx_IN_current = inpolygon(PClast(:,1),PClast(:,2),CAfit_rad{irad}(:,1),CAfit_rad{irad}(:,2));
        if irad ~= 1
            idx_IN_allprevious = idx_IN_allprevious|idx_IN_rad{irad-1};
        end
        idx_IN_rad{irad} = idx_IN_current & ~idx_IN_allprevious; %& idx_frontal_last;
    end
    CAfit =[CAfit,repelem(mean(CAlast(:,3),'omitnan'),length(CAfit),1)];
    PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
    ca_idx_current = T{mm}.mask_in_border(:,frame_target);
    % graph
    
    for irad = 1:Nrad_segmentation
        fig_name = sprintf('map annulus radial num %d',irad); 
        fig = newfig(fig_name); Nrow = 1; Ncol = 1; 
        set(fig,'Position',[17.0000   14.0000    4.8281    4.1240]); 
        h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
        ax_idx = 1; 
        scatter3(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),PCcurrent(:,3),2,[1 1 1].*0.5,"filled");
        scatter3(h_ax(ax_idx),PCcurrent(idx_IN_rad{irad},1),PCcurrent(idx_IN_rad{irad},2),PCcurrent(idx_IN_rad{irad},3),10,'k',"filled");
        plot3(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),CAfit(:,3),'-','Color',[0.8500 0.3250 0.0980])
        set(h_ax(ax_idx),'XLim',[-15,15],'YLim',[-10,15],...
            'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','View', [0 -90],'visible','off');
        colormap(h_ax(ax_idx),'jet');
    end

end


%% Visualisation / interpretation of results. 
list_face_target_num = [1,2]; 
list_rad_num = 1:Nrad_segmentation; 
list_force_target = [0.01,0.02,0.03,0.04,0.05,0.10,0.20,0.50,1:1:5];

magnification = 5; %magnification / amplification for visual purpose
square_w = 1; 

fig_name = sprintf('Visu schematic deformation %s',SUBJECT); 
fig = newfig(fig_name); 
set(fig,'Position',[1.3494   14.3140    8.4931   11.0596]); 
Nrow = 1; Ncol = 1; 
h_ax = subplot_ax(Nrow,Ncol);
% for mm = 1%:Ntrial
%     itrial = itrial + 1; 
%     ax_idx = 1;  
%     FCmat_cur = FCmat_cum{mm};
%     for kk = 1:length(list_frame_num)
%         FCmat_epc1 = FCmat_cur(list_frame_num(kk),list_rad_num+0*Nrad_segmentation);
%         FCmat_epc2 = FCmat_cur(list_frame_num(kk),list_rad_num+1*Nrad_segmentation);
%         
%         a = square_w; 
%         b = square_w; 
%         h = square_w; 
%         h_cum = 0; 
%         ii = 0; 
%         for irad = list_rad_num
%             ii = ii + 1;  
%             h = square_w*(1+FCmat_epc1(ii)/1e2*magnification);
%             b = square_w*(1+FCmat_epc2(ii)/1e2*magnification);
%             XYtrap = trapezium(a,b,h);
%             patch(h_ax(ax_idx),XYtrap(:,1)+kk*2,XYtrap(:,2)+h_cum,'r','FaceAlpha',0.1);
%             a = b; 
%             h_cum = h_cum + h; 
%         end
%         xline(h_ax(ax_idx),kk*2+square_w);
%     end
% end
list_rad_num = 1:Nrad_segmentation;
for mm = 1:Ntrial
    nf_format = T{mm}.nf;
    nf_format = nf_format(1:find(nf_format>NFCOND_SET,1,'first'));
    list_frame_num = find_closest_value(nf_format,list_force_target);
    frame_target = list_frame_num(1); 
    itrial = itrial + 1; 
    ax_idx = 1;  
    FCmat_cur = FCmat_cum{mm};
    for kk = 1:length(list_frame_num)
        FCmat_epc1 = FCmat_cur(list_frame_num(kk),list_rad_num+0*Nrad_segmentation);
        FCmat_epc2 = FCmat_cur(list_frame_num(kk),list_rad_num+1*Nrad_segmentation);
        
        
        ii = ii + 1; h_cum = 0;  
        p1 = mean(FCmat_epc1(1:Nrad_in));
        p2 = mean(FCmat_epc2(1:Nrad_in));
        
        if kk == length(list_frame_num)
            fa = 0.2;
        else
            fa = 0.0; 
        end
        a = square_w*(1+p2/1e2*magnification); 
        h = square_w*(1+p1/1e2*magnification);
        b = a;
        XYtrap = trapezium(a,b,h);
        patch(h_ax(ax_idx),XYtrap(:,1),XYtrap(:,2)+h_cum,'k','FaceAlpha',fa);
        ii = ii + 1; h_cum = 0;%h_cum+h;  
        p1 = mean(FCmat_epc1(Nrad_in+1:Nrad_segmentation));
        p2 = mean(FCmat_epc2(Nrad_in+1:Nrad_segmentation));

        a = square_w*(1+p2/1e2*magnification); 
        h = square_w*(1+p1/1e2*magnification);
        b = a;
        XYtrap = trapezium(a,b,h);
        patch(h_ax(ax_idx),XYtrap(:,1)+square_w*2,XYtrap(:,2)+h_cum,'k','FaceAlpha',fa);
    end
    text(h_ax(ax_idx),0,0,...
        'Inside','VerticalAlignment','top',...
        'HorizontalAlignment','center',...
        'FontSize',15);
    
    text(h_ax(ax_idx),square_w*2,0,...
        'Outside','VerticalAlignment','top',...
        'HorizontalAlignment','center',...
        'FontSize',15);
end
set(h_ax(ax_idx),'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','visible','off');


list_rad_num = 1:Nrad_segmentation; 
list_force_target = 0:0.1:5;
color_frame = jet(length(list_force_target));%(1-linspace(0,1,length(list_frame_num)))'.*[1 1 1];% 
color_finger_part_segmentation = [0.8500 0.3250 0.0980]; %orange

fig_name = sprintf('Visu schematic deformation 2 %s x%d',SUBJECT,magnification); 
fig = newfig(fig_name); 
set(fig,'Position',[1.5346   18.6796   48.5775    7.5671]); 
Nrow = Ntrial; Ncol = length(list_rad_num); 
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
itrial = 0; 
for mm = 1:Ntrial
    itrial = itrial + 1; 
    nf_format = T{mm}.nf;
    nf_format = nf_format(1:find(nf_format>NFCOND_SET,1,'first'));
    list_frame_num = find_closest_value(nf_format,list_force_target);
    frame_target = list_frame_num(1); 
    FCmat_cur = FCmat_cum{mm};
    for irad = 1:length(list_rad_num)
        ax_idx = irad+Ncol*(itrial-1);  
        for kk = 1:length(list_frame_num)
            FCmat_epc1 = FCmat_cur(list_frame_num(kk),list_rad_num(irad)+0*Nrad_segmentation);
            FCmat_epc2 = FCmat_cur(list_frame_num(kk),list_rad_num(irad)+1*Nrad_segmentation);
            p1 = FCmat_epc1;
            p2 = FCmat_epc2;

            a = square_w*(1+p2/1e2*magnification); 
            h = square_w*(1+p1/1e2*magnification);
            b = a;%square_w*(1+p2/1e2*magnification);
            XYtrap = trapezium(a,b,h);
%             if kk == length(list_frame_num)
                patch(h_ax(ax_idx),XYtrap(:,1),XYtrap(:,2),nf_format(list_frame_num(kk)),'FaceAlpha',1,'EdgeColor','none');
%             else
%                 patch(h_ax(ax_idx),XYtrap(:,1),XYtrap(:,2),nf_format(list_frame_num(kk)),'FaceAlpha',0.0);
%             end
        end
        if irad == 1 
            rad_combined_number = 0.1; L = square_w+2*rad_combined_number; 
            p_origin = [-rad_combined_number,-rad_combined_number];
            quiver(h_ax(ax_idx),p_origin(1),p_origin(2),L,0,'Color','k','ShowArrowHead','on');
            quiver(h_ax(ax_idx),p_origin(1),p_origin(2),0,L,'Color','k','ShowArrowHead','on');
            text(h_ax(ax_idx),p_origin(1)+L,p_origin(2),...
                sprintf('circum'),'VerticalAlignment','top',...
                'HorizontalAlignment','right',...
                'FontSize',15);
            text(h_ax(ax_idx),p_origin(1),p_origin(2)+L,...
                sprintf('radial'),'VerticalAlignment','top',...
                'HorizontalAlignment','right',...
                'FontSize',15);
        end
        text(h_ax(ax_idx),0,square_w+rad_combined_number,...
            sprintf('#%d',irad),'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'FontSize',15);
        set(h_ax(ax_idx),'XLim',[-square_w/2,square_w/2],'YLim',[0,square_w],...
            'clim',[0,5]);
        colormap(h_ax(ax_idx),flipud(jet));
    end    
    text(h_ax(ax_idx),square_w/2,0,...
        sprintf('Magn : x%d',magnification),'FontSize',20,...
        'VerticalAlignment','top','HorizontalAlignment','right');
    
%     ax_idx = Ncol;
%     colormap(h_ax(ax_idx),'jet');
%     h_cb = colorbar(h_ax(ax_idx),...
%         'Ticks',list_force_target,...
%         'TickLabels',num2str(list_force_target'),...
%         'FontWeight','normal',...
%         'FontSize',15,...
%         'FontName','Arial',...
%         'Location','east',...
%         'AxisLocation','out'); 
%     h_cb.Title.String=sprintf('Force (N)'); 
%     set(h_ax(ax_idx),'clim',[0,list_force_target(end)],'Color','none'); 
    
    % map 
%     if itrial == 1
%         PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
%         PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
%         ca_idx_current = T{mm}.mask_in_border(:,frame_target);
% 
%         % graph 
%         for irad = 1:length(list_rad_num)
%             ax_idx = 1*Ncol+irad;
%             scatter(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),1,[1 1 1].*0.8,"filled");
%             scatter(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),2,color_finger_part_segmentation,"filled");
%     %         plot(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
%             set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
%                 'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
%             colormap(h_ax(ax_idx),'jet');
%         end
%         qw = cell(1,1); 
%         qw{1} = fill(h_ax(ax_idx),[NaN NaN],[NaN NaN],color_finger_part_segmentation,'EdgeColor','none');
%         lgd_label = {"Measured part"};
%         legend(h_ax(ax_idx),[qw{:}],lgd_label,'Location','southeast','box','off')
%     end
end

set(h_ax,'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','visible','off');



%% CF search 
face_target = 1:2; 
trial_target = [1,4]; 

mvaluecf1 = NaN(length(trial_target),1);
mvaluestrain = NaN(length(trial_target),length(face_target));
title_fig = 'Search for friction relation with loading 1'; 
fig = newfig(title_fig);
set(fig,'Unit','Normalized','Position',[0.3346    0.0287    0.1673    0.8926]);
Nrow = 5; Ncol = 1; merge_element = {}; 
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on'); grid(h_ax,'on');

for iface = face_target 
    itrial = 0;   
    for mm = trial_target
        itrial = itrial + 1; 
        nf_format = T{mm}.nf;
        tf_format = T{mm}.tf;

        idx_final_force = find(nf_format>NFCOND_SET,1,'first'); 
        nf_format = nf_format(1:idx_final_force);
        tf_format = tf_format(1:idx_final_force);
%         idx_final_force = find(nf_raw(:,mm)>1,1,'first'); 
%         nf_format = nf_raw(1:idx_final_force,mm);
%         tf_format = tf_raw(1:idx_final_force,mm);
        if iface == 1
            ax_idx = 1; 
            plot(h_ax(ax_idx),nf_format,nf_format,'.-','Color',C(mm,:),'MarkerSize',10);
            ylabel(h_ax(ax_idx),'NF (N)');
            
            ax_idx = 2; 
            plot(h_ax(ax_idx),nf_format,tf_format,'.-','Color',C(mm,:),'MarkerSize',10);
            ylabel(h_ax(ax_idx),'TF (N)');

            ax_idx = 3; 
            ratio_force = tf_format./nf_format; 
            hp = plot(h_ax(ax_idx),nf_format,ratio_force,'.-','Color',C(mm,:),'MarkerSize',10);
            set(h_ax(ax_idx),'YLim',[0,2]);
            xlabel(h_ax(ax_idx),'Force (N)');
            ylabel(h_ax(ax_idx),'TF/NF (-)');
            
            idx_start_force = find(nf_format>0.05,1,'first');
            idx_1N_force = find(nf_format>1,1,'first');
            temp = ratio_force; mask = false(length(ratio_force),1); mask(idx_start_force:idx_1N_force)=true; 
            temp(~mask) = NaN; 
            [mvaluecf1(itrial),idx] = max(temp);
%             mvaluecf1(mm) = temp(idx_1N_force);
%             idx = idx_1N_force;
            
            plot(h_ax(ax_idx),nf_format(idx),ratio_force(idx),'o',...
                'Color',C(mm,:),'MarkerFaceColor',C(mm,:),'MarkerEdgeColor','k','MarkerSize',8);
        end

        seg_target = 3:Nrad_in;
        FCmat_cur = FCmat_cum{mm};
        idx_1N_force = find(nf_format>NFCOND_SET*0.95,1,'first');
        mvaluestrain(itrial,iface) = mean(FCmat_cur(idx_1N_force,seg_target+(iface-1)*Nrad_segmentation));
        ax_idx = 3+iface; 
        plot(h_ax(ax_idx),mvaluecf1(itrial),mvaluestrain(itrial,iface),'.','Color',C(mm,:));
%         plot(h_ax(ax_idx),mvaluecf2(mm,iface),mvaluestrain(mm,iface),'.','Color',C(mm,:));
    end
    ax_idx = 3+iface;
    plot(h_ax(ax_idx),mvaluecf1,mvaluestrain(:,iface),'k-');
%     plot(h_ax(ax_idx),mvaluecf2(:,iface),mvaluestrain(:,iface),'k-');
    xlabel(h_ax(ax_idx),'max(TF/NF(N))');
    ylabel(h_ax(ax_idx),FACEMEASURE{iface});
end

%     ylim(h_ax(ax_idx),[-10,-5]); 
legend(h_ax(1),"trial #"+num2str((TRIAL_TARGET)')); 
set(h_ax,'FontSize',12,'clipping','on','Color','None');
% xlim(h_ax(ax_idx),[0.5,1]); 

% Find a correlation between coefficient of friction and 
%%
figure; 
h_ax = subplot_ax(2,1); hold(h_ax,'on'); 
for mm = 1:Ntrial
    idx_start_force = find(nf_raw(:,mm)>0.05,1,'first');
    ax_idx = 1; 
    plot(h_ax(ax_idx),nf_raw(idx_start_force:end,mm));
    plot(h_ax(ax_idx),tf_raw(idx_start_force:end,mm));
    ax_idx = 2; 
    plot(h_ax(ax_idx),tf_raw(idx_start_force:end,mm)./nf_raw(idx_start_force:end,mm));
    ylim(h_ax(ax_idx),[0,4]) 
end
% title_fig = 'Search for friction relation with tangential 2'; 
% fig = newfig(title_fig);
% Nrow = 2; Ncol = 1; merge_element = {}; 
% h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on'); 
% 
% itrial = 0; 
% for mm = 1:Ntrial
%     itrial = itrial + 1; 
%     
%     nf_format = T{mm}.nf;
%     tf_format = T{mm}.tf;
%     
%     idx_tang_loading = 151; 
%     nf_format = nf_format(idx_tang_loading:end);
%     tf_format = tf_format(idx_tang_loading:end);
%     
%     ax_idx = 1; 
%     plot(h_ax(ax_idx),1:length(nf_format),tf_format,'.-');
%     
%     ax_idx = 2; 
%     plot(h_ax(ax_idx),1:length(nf_format),tf_format./nf_format,'.-');
%     set(h_ax(ax_idx),'YLim',[-3,3]);
% end
% legend(h_ax(1),num2str((1:Ntrial)')); 
% set(h_ax,'FontSize',12,'clipping','off','Color','None');
    



%% Maximum rate for multiple rad with force parameter

list_force_target = [0,0.02,0.03,0.04,0.05,0.10,0.20,0.50,1:1:5];
% list_frame_num = 27:1:90;
ylim_set = [-50,50]; 
face_target = 1; Nface_target = length(face_target); 
facemaxrate = zeros(Ntrial,Nrad_segmentation,Nface_target); 
idx_max_rate = zeros(Ntrial,Nrad_segmentation,Nface_target); 
color_face = [[1 1 1].*0.7;...
              [1 1 1].*0.0];

fig = newfig(sprintf('Multiple segmentation : Max rate vs spatial position')); 
set(fig,'Position',[3.7306   10.1865   13.4144   13.9435]); 
Nrow = 7; Ncol = Nrad_segmentation; 
merge_element = {1:5*Ncol};
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');
itrial = 0; 
ph = cell(1,2); 
itrial = 0; 
for mm = 1:Ntrial
    itrial = itrial+1; 
    FCmat_cur = FCmat_rate{mm}; 
    
    nf_format = T{mm}.nf;
    nf_format(find(nf_format>5,1,'first')+1:end) = NaN;
    list_frame_num = find_closest_value(nf_format,list_force_target)';
    for iface = face_target
        for irad = 1:Nrad_segmentation 
            FCmat_current = FCmat_cur(list_frame_num,irad+(iface-1)*Nrad_segmentation);
            [~,idx] = max(abs(FCmat_current)); 
            idx_max_rate(mm,irad,iface) = idx; 
            facemaxrate(mm,irad,iface) = FCmat_current(idx);
        end
    end
end

%graph
ax_idx = 1; 
yline(h_ax(ax_idx),0,'Color',[1 1 1].*0.5); 
plot(h_ax(ax_idx),repelem(1:Nrad_segmentation,2,1),[ylim_set(1).*ones(1,Nrad_segmentation);...
    satur(min(FCmat_cur(list_frame_num,1:Nrad_segmentation),[],1),'level',ylim_set(1),'method','low')],'-','Color',[1 1 1].*0.5);

for iface = face_target
    for irad = 1:Nrad_segmentation 
        scatter(h_ax(ax_idx),irad,FCmat_current(idx),30,nf_all(list_frame_num(idx),mm),'filled');

%         text(h_ax(ax_idx),irad,facemaxrate(irad,iface),...
%             sprintf('%1.2f',nf_all(list_frame_num(idx_max_rate(irad,iface)),mm)),...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','top',...
%             'FontSize',8);
    end
    size(1:Nrad_segmentation)
    size(facemaxrate(:,iface))
    ph{iface}=bar(h_ax(ax_idx),1:Nrad_segmentation,facemaxrate(1,:,iface),...
        'FaceColor',color_face(iface,:),'EdgeColor','None'); 
end
uistack([ph{2}],'bottom');
uistack([ph{1}],'bottom');

colormap(h_ax(ax_idx),'jet'); %[1 1 1].*color_scale_grey'); 
h_cb = colorbar(h_ax(ax_idx)); set(h_cb,'Location','eastoutside');
h_cb.Title.String=sprintf('Force (N)'); 
set(h_ax(ax_idx),'clim',[0,5],'Color','none'); 
h_ax(ax_idx).YLabel.String = 'Rate (%)'; 
set(h_ax(ax_idx),'XLim',[0,Nrad_segmentation+1],'YLim',ylim_set,'Clipping','off');

if itrial == 1
    % map 
    PCcurrent = T{mm}.FC.xyz(:,:,frame_target);
    PCcurrent_segmented = squeeze(FCmat_xyz_all{mm}(:,frame_target,:,:));
    ca_idx_current = T{mm}.mask_in_border(:,frame_target);

    % graph 
    for irad = 1:Nrad_segmentation
        ax_idx = 5*Ncol+irad;
        if mod(ax_idx,2)==0 
            ax_idx = ax_idx+Ncol;
        end
        scatter(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),1,[1 1 1].*0.8,"filled");
        scatter(h_ax(ax_idx),PCcurrent_segmented(:,1,irad),PCcurrent_segmented(:,2,irad),2,[0.8500 0.3250 0.0980],"filled");
%         plot(h_ax(ax_idx),CAfit(:,1),CAfit(:,2),'-','Color',[1 1 1].*0.5,'LineWidth',1)
        set(h_ax(ax_idx),'XLim',[-7,7],'YLim',[0,5],...
            'FontSize',12,'clipping','off','DataAspectRatio',[1 1 1],'Color','None','view',[0 -90],'visible','off');
        colormap(h_ax(ax_idx),'jet');
    end

    set(h_ax,'visible','off');
    set(h_ax(1),'visible','on');
    set(get(h_ax(1), 'XAxis'), 'Visible', 'off');
    set(h_ax,'FontSize',12); 
end


%% Maximum rate for multiple rad with force parameter

list_frame_num = 27:1:95;
ylim_set = [-50,50]; 
facemaxrate = zeros(length(list_frame_num),2); 
idx_max_rate = zeros(length(list_frame_num),2); 
color_face = [[1 1 1].*0.7;...
              [1 1 1].*0.0];

fig = newfig(sprintf('Multiple segmentation : Max rate vs force')); 
set(fig,'Position',[3.7306   10.1865   13.4144   13.9435]); 
Nrow = 1; Ncol = 1; 
merge_element = {};
h_ax = subplot_ax(Nrow,Ncol,'merge',merge_element); hold(h_ax,'on');
itrial = 0; 
ph = cell(1,2); 
for mm = 1:Ntrial
    itrial = itrial+1; 
    ax_idx = 1; 
    FCmat_cur = FCmat_rate{mm}; 
    
    yline(h_ax(ax_idx),0,'Color',[1 1 1].*0.5); 
    
    for iface = 1:2
        frame_count = 0; 
        for iframe = list_frame_num 
            frame_count = frame_count+1; 
            FCmat_current = FCmat_cur(iframe,(1:Nrad_segmentation)+(iface-1)*Nrad_segmentation);
            [~,idx] = max(abs(FCmat_current)); 
            idx_max_rate(frame_count,iface) = idx; 
            facemaxrate(frame_count,iface) = FCmat_current(idx);
        end
        nf_current = nf_all(list_frame_num,mm);
        ph{iface}=plot(h_ax(ax_idx),nf_current,facemaxrate(:,iface),...
            'Color',color_face(iface,:)); 
    end
    uistack([ph{2}],'bottom');
    uistack([ph{1}],'bottom');
    
    set(h_ax(ax_idx),'XLim',[0,5],'YLim',ylim_set,'Clipping','off');

    
end





