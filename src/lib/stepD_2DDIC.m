
function varargout = stepD_2DDIC(varargin)
%% Parse input arguments
p = inputParser;
p.addParameter('baseDataPath',[]);
p.addParameter('baseResultPath',[]);
p.addParameter('subject',[]);
p.addParameter('material',[]);
p.addParameter('trial',[]);
p.addParameter('stereopair',[]);
p.addParameter('phase',[]);
p.addParameter('jump',[]);
p.addParameter('idxstart_set',[]);
p.addParameter('idxend_set',[]);
p.addParameter('showvisu',0);
p.addParameter('savedata',1);
p.addParameter('reftrial_setmanual',[]);
p.addParameter('param_filt_im',[25 300]);
p.addParameter('automatic_process', true);

p.parse(varargin{:});
baseDataPath = p.Results.baseDataPath;
baseResultPath = p.Results.baseResultPath;
subject = p.Results.subject;
material = p.Results.material;
trial = p.Results.trial;
stereopair = p.Results.stereopair;
phase = p.Results.phase;
idxstart_set = p.Results.idxstart_set;
idxend_set = p.Results.idxend_set;
jump = p.Results.jump;
showvisu = p.Results.showvisu;
savedata = p.Results.savedata;
reftrial_setmanual = p.Results.reftrial_setmanual;
param_filt_im = p.Results.param_filt_im;
automatic_process = p.Results.automatic_process;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TEST PARAMETERS : TO BE COMMENTED !!!!
% baseDataPath = pwd;
% subject = 'S02';
% material = 'glass';
% trial = '003';
% stereopair = 1;
% phase = 'slide1';
% jump = 5;
% showvisu = 1;
% savedata = 1;
% idxstart_set=[];
% idxend_set=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("-------------------------------------------\n"); 
fprintf("-------------------------------------------\n"); 
fprintf("Digital Image Correlation analysis launch\n"); 

%Global constant 
TRUE_FPS = 50; 
if str2double(cell2mat(regexp(subject, '\d+', 'match')))<8
    LIMIT_GRAYSCALE = 70;%70;
else
    LIMIT_GRAYSCALE = 100;
end

% Filtering parameters
im_filter_mode = true; 

% Ncorr DIC parameters
analysis_direction = 'regular'; 
subset_radius_ncorr_tracking = 40; %40%careful :35
subset_radius_ncorr_matching = 60; %50%S02 : 60 
subset_spacing = 10;%10 %20
cutoff_tracking = 1e-5; 
cutoff_matching = 1e-5; 
number_iteration_solver = 100; 
number_threads = 1; 
high_strain_analysis = 1; 
seed_propagation = 'seed'; 
auto_ref_change = 1; 
step_ref_change = 10;

% ensure output path exists
outputPath = fullfile(baseResultPath,subject,material,trial,phase); %,sprintf('filt_%d_%d',param_filt_im(1),param_filt_im(2)));
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end

% load protocol
protocolPath = fullfile(fullfile(baseDataPath,"rawdata",subject,"speckles",...
    material,"protocol",sprintf('*.mat')));
S = dir(protocolPath);
if isempty(S)
    error('Error: Protocol not found.');
end
p = load(fullfile(S.folder,S.name)); protocol = p.cond;

dircond = protocol.table(:,strcmp(protocol.titles,'dir'));
nfcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'nf')));
spdcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'spd')));
repcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'rep')));
spddxlcond = repelem([0.08],length(repcond)); 


% reference trial for ROI tracing
if (strcmp(phase,"loading") || nfcond(str2double(trial)) == 1)
    %REF for all when phase loading
    %REF for 1N when in slide1 phase
    reftrial = sprintf("%03d",find(nfcond == 1&...
        strcmp(dircond,"Ubnf"),1)); 
elseif strcmp(phase,"slide1") && nfcond(str2double(trial)) == 5 
    %REF for 5N when in slide1 phase
    reftrial = sprintf("%03d",find(nfcond == 5 &...
        strcmp(dircond,"Ubnf"),1));
else
    fprintf("error : no ref trial assigned\n"); 
    return; 
end

%manual set
if ~isempty(reftrial_setmanual)
    reftrial = reftrial_setmanual;
end

%pair order for 3D reconstruction 
if strcmp(dircond(str2double(trial)),"Ubnf")
    pairOrder = [2,1];
    pairForced = true;
elseif strcmp(dircond(str2double(trial)),"Rbnf") 
    pairOrder = [1,2];
    pairForced = true;
else 
    pairOrder = [1,2];
    pairForced = false; 
end 

%% parameter structure
%filename
roifile = fullfile(baseResultPath,subject,material,...
    "REF_MASK_"+reftrial+"_"+phase+"_pair"+num2str(stereopair)+".mat"); 
matchingfile = fullfile(outputPath,...
    "MATCHING2"+reftrial+"_pair"+num2str(stereopair)+".mat");
seedfile = fullfile(baseResultPath,subject,material,...
    "REF_SEED_"+reftrial+"_"+phase+"_pair"+num2str(stereopair)+".mat"); 

base_parameters = struct(...
    'baseDataPath',baseDataPath,...
    'baseResultPath',baseResultPath,...
    'subject',subject,...
    'material',material,...
    'trial',trial,...
    'stereopair',stereopair,...
    'phase',phase,...
    'jump',jump,...
    'outputPath',outputPath,...
    'reftrial',reftrial,...
    'roifile',roifile,...
    'matchingfile',matchingfile,...
    'seedfile',seedfile,...
    'idxstart_set',idxstart_set,...
    'idxend_set',idxend_set,...
    'limit_grayscale',LIMIT_GRAYSCALE); 

step1_parameters = struct( ...
    'type',analysis_direction,...
    'radius',subset_radius_ncorr_tracking,...
    'spacing',subset_spacing, ...
    'cutoff_diffnorm',cutoff_tracking, ...
    'cutoff_iteration',number_iteration_solver, ...
    'total_threads',number_threads, ...
    'stepanalysis_params', struct( ...
                            'enabled',high_strain_analysis, ...
                            'type',seed_propagation, ...
                            'auto',auto_ref_change, ...
                            'step',step_ref_change) ...
    );
step1_2_parameters = struct( ...
    'type',analysis_direction,...
    'radius',subset_radius_ncorr_matching,...
    'spacing',subset_spacing, ...
    'cutoff_diffnorm',cutoff_matching, ...
    'cutoff_iteration',number_iteration_solver, ...
    'total_threads',number_threads, ...
    'stepanalysis_params', struct( ...
                            'enabled',high_strain_analysis, ...
                            'type',seed_propagation, ...
                            'auto',auto_ref_change, ...
                            'step',step_ref_change) ...
    );

step2_parameters = step1_parameters;
fprintf('global parameters set\n');


%% Read data
% read images from video file and retreive camera number according to
% stereopair given.
tic; 
[cam_first_raw, cam_second_raw, cam_1, cam_2] = import_vid(baseDataPath,...
                                       'subject',subject,...
                                       'material',material,...
                                       'trial',trial,...
                                       'stereopair',stereopair,...
                                       'phase',phase,...
                                       'idxstart_set',idxstart_set,...
                                       'idxend_set',idxend_set,...
                                       'framejump',jump);
elapsedTime = toc; 
fprintf("reading done in %1.1fs\n",elapsedTime); 

%% filtering images : use all images
%limit_grayscale_level = 75; 
if strcmp(phase,"slide1")
    cam_first_raw = cam_first_raw(:,:,1:end/2+5);
    cam_second_raw = cam_second_raw(:,:,1:end/2+5);
end

%saturation of images
cam_first_satur = satur(cam_first_raw,'level',LIMIT_GRAYSCALE); 
cam_second_satur = satur(cam_second_raw,'level',LIMIT_GRAYSCALE);


%% Initialization of ROI - SEED - MATCHING PAIR 
% Load ROI
if ~exist(roifile,'file')    
    draw_ref_roi(base_parameters);
end

% Matching is done between a reference trial to the current one to skip
% making the roi for each cameras at each trial. It is done over the
% saturated images 
% verify that the matching exist and if not, ncorr is run to acquire 
% displacement field data from 
if ~exist(matchingfile,'file')    
    ncorr_matching2ref(cam_first_satur(:,:,1),base_parameters,step1_2_parameters);
end

% compute reformatted roi
matching = load(matchingfile); 
refmask_REF = matching.reference_save(1).roi.mask; 
refmask_trial = matching.current_save(1).roi.mask; 
disp('--> STEP : ROI loaded and formatted'); 

% SEED placement
if ~exist(seedfile,'file')
    %manual seed placement in the reference trial 
    draw_ref_seed(base_parameters,'roi',refmask_REF);
end

%load seed position
seed_point = load(seedfile);

%mapping coordinate
ref_seed_point.pw = seed_point.seed_point; %pixel world
ref_seed_point.sw = map_pixel2subset(ref_seed_point.pw,subset_spacing); %subset world

% mapping : format seed point for ncorr1 and ncorr12
U_mapped = matching.data_dic_save.displacements(1).plot_u_ref_formatted/(subset_spacing+1); 
V_mapped = matching.data_dic_save.displacements(1).plot_v_ref_formatted/(subset_spacing+1); 

initial_seed_point_set1.sw = map_pointcoordinate(ref_seed_point.sw,{U_mapped,V_mapped}); 
initial_seed_point_set1.pw = map_subset2pixel(initial_seed_point_set1.sw,subset_spacing);

disp('--> STEP : SEED loaded and formatted');
%%
if im_filter_mode == true 
    % We reuse the mask that has been processed to keep bounds of gs
    % intensity over the roi of the first cam at initial step. 
    [cam_first,gsboundaries] = filter_like_ben(cam_first_satur,'mask',refmask_trial,'paramfilt',param_filt_im);
    % Keep the same bounds (should be a good approximation) 
    cam_second = filter_like_ben(cam_second_satur,'gsbound',gsboundaries,'paramfilt',param_filt_im); 
    disp('--> STEP : filtering done');
else 
    cam_first = cam_first_raw; 
    cam_second = cam_second_raw; 
    disp('--> STEP : raw data used');

end
%% trial info
idxframe = idxstart_set:jump:idxend_set;
actual_fps_meas = TRUE_FPS/jump; 
fprintf("Trial information summary \'%s\' %s %s %s\n   > Dir: %s\n   > Force: %dN\n   > Spd: %dmm/s\n   > Rep: %d\n   > REFtrial: \'%s\'\n   > NbrFr: %d\n   > Frame : %s\n   > Actual FPS : %d\n",...
    trial,...
    phase,...
    material,...
    subject,...
    dircond{str2double(trial)},...
    nfcond(str2double(trial)),...
    spdcond(str2double(trial)),...
    repcond(str2double(trial)),...
    reftrial,...
    size(cam_first,3),...
    mat2str(idxstart_set:jump:idxend_set),...
    actual_fps_meas);

savefileName = fullfile(outputPath,sprintf('dic_info_data_target_pair%d.mat',stereopair));
save(savefileName,'actual_fps_meas','idxframe'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ncorr analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATCHING STEP - RUN #1
tic; 
step1_2_parameters.initial_seed = initial_seed_point_set1.pw; 
% We keep the ROI of the ncorr matching step

%
siz = size(cam_first); 
input1 = cam_first_satur(:,:,1); 
input2 = zeros(siz(1),siz(2),2,'uint8'); 
input2(:,:,1) = cam_second_satur(:,:,1); 
input2(:,:,2) = cam_first_satur(:,:,1); %needed to avoid a bug in ncorr
[h12,file_logic] = ncorr_dic_rewrited('cam_number',[cam_1,cam_2],...
    'cam_data_ref',input1,...
    'cam_data_cur',input2,...
    'mask',refmask_trial,...
    'automatic_process',automatic_process,...
    'step_param',step1_2_parameters,...
    'base_param',base_parameters);

% Retreive info for next step 
refmask_trial_matched = h12.current(1).roi.mask;
% mapping : format seed point for ncorr2
U_mapped = h12.data_dic.displacements(1).plot_u_ref_formatted/(subset_spacing+1); 
V_mapped = h12.data_dic.displacements(1).plot_v_ref_formatted/(subset_spacing+1); 

initial_seed_point_set2.sw = map_pointcoordinate(initial_seed_point_set1.sw,{U_mapped,V_mapped}); 
initial_seed_point_set2.pw = map_subset2pixel(initial_seed_point_set2.sw,subset_spacing);

if ~file_logic
    close(h12.handles_gui.figure); 
end
clear('h12');
elapsedTime = toc; 
fprintf("--> STEP : Ncorr matching 1-2 done in %1.1fs\n",elapsedTime);


%% Visualisation 
if showvisu
    siz = size(cam_first);
    cam_first_1 = cam_first(:,:,1); 
    cam_second_1 = cam_second(:,:,1); 
    cam_first_1(~refmask_trial) = 0.25*cam_first_1(~refmask_trial);
    cam_second_1(~refmask_trial_matched) = 0.25*cam_second_1(~refmask_trial_matched);
    h1 = player([repmat([reshape(cam_first_1,[siz(1),siz(2)]),reshape(cam_second_1,[siz(1),siz(2)])],1,1,siz(3));...
                cam_first,cam_second],'fps',TRUE_FPS/jump,'movie',0); hold on; 
    drawpoint('Position',[initial_seed_point_set1.pw(1),initial_seed_point_set1.pw(2)]); hold on; 
    drawpoint('Position',[initial_seed_point_set2.pw(1)+siz(2),initial_seed_point_set2.pw(2)]); 
      
    pause(30); 
    
    % if not already closed, close player
    isDeletedObj = isobject(h1) & ~isgraphics(h1);
    if ~isDeletedObj
         close(h1); 
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACKING STEP 1

tic; 
step1_parameters.initial_seed = initial_seed_point_set1.pw; 
% We keep the ROI of the ncorr matching step
[h1,file_logic] = ncorr_dic_rewrited('cam_number',cam_1,...
    'cam_data_ref',cam_first(:,:,1),...
    'cam_data_cur',cam_first(:,:,1:end),...
    'mask',refmask_trial,...
    'automatic_process',automatic_process,...
    'step_param',step1_parameters,...
    'base_param',base_parameters);

if ~file_logic
    close(h1.handles_gui.figure); 
end
clear('h1');
elapsedTime = toc; 
fprintf("--> STEP : Ncorr 1 done in %1.1fs\n",elapsedTime); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACKING STEP 2

tic; 
step2_parameters.initial_seed = initial_seed_point_set2.pw;
% We keep the ROI of the ncorr matching step
[h2,file_logic] = ncorr_dic_rewrited('cam_number',cam_2,...
    'cam_data_ref',cam_second(:,:,1),...
    'cam_data_cur',cam_second(:,:,2:end),...
    'mask',refmask_trial_matched,...
    'automatic_process',automatic_process,...
    'step_param',step2_parameters,...
    'base_param',base_parameters);

if ~file_logic
    close(h2.handles_gui.figure); 
end
clear('h2');
elapsedTime = toc; 
fprintf("--> STEP : Ncorr 2 done in %1.1fs\n",elapsedTime); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format output

step2_dic_finish(outputPath,...
                'cam_1',cam_1,...
                'cam_2',cam_2);
if savedata == 0
    delete(fullfile(outputPath,'ncorr1.mat'));
    delete(fullfile(outputPath,'ncorr12.mat'));
    delete(fullfile(outputPath,'ncorr2.mat'));
end

disp('--> STEP : Ncorr analysis completed');

if nargout >= 1
    varargout{1} = outputPath; 
end
if nargout >= 2
    varargout{2} = pairOrder; 
end
if nargout >= 3
    varargout{3} = pairForced; 
end

%% Visu 
% if showvisu == 1
%     bool_point = true; 
%     factor = 4; 
%     while bool_point
%         siz = size(cam_first_satur); 
%         fig = figure; imshow([cam_first_satur(:,:,1),cam_second_satur(:,:,1)],[]);
%         [X,Y] = ginput(2);
%         close(fig);
%         c1 = round([X(1),Y(1)]); c2 = round([X(2)-siz(2),Y(2)]);
%         R = subset_radius_ncorr_tracking; 
% 
%         [x1,y1] = circle_coordinates(c1,R); 
%         [xy1] = findboundary([x1,y1],'connex_factor',1);
%         x1 = xy1(:,1); y1 = xy1(:,2); 
%         mask1 = poly2mask(x1,y1,siz(1),siz(2)); 
%         maxxy=max(xy1); minxy=min(xy1);
%         mask1_small = mask1(minxy(2):maxxy(2),minxy(1):maxxy(1));
% 
%         [x2,y2] = circle_coordinates(c2,R); 
%         [xy2] = findboundary([x2,y2],'connex_factor',1); 
%         x2 = xy2(:,1); y2 = xy2(:,2); 
%         mask2 = poly2mask(x2,y2,siz(1),siz(2));
%         maxxy=max(xy2); minxy=min(xy2);
%         mask2_small = mask2(minxy(2):maxxy(2),minxy(1):maxxy(1)); 
% 
%         im_first_raw = cam_first_raw(:,:,1); 
%         im_second_raw = cam_second_raw(:,:,1); 
%         im_first = cam_first_satur(:,:,1); 
%         im_second = cam_second_satur(:,:,1);
% 
%         %highlight circle
%         im_first_raw(~mask1) = im_first_raw(~mask1)/factor;
%         im_second_raw(~mask2) = im_second_raw(~mask2)/factor;
%         im_first(~mask1) = im_first(~mask1)/factor;
%         im_second(~mask2) = im_second(~mask2)/factor;
% 
%         %zoom
%         z = zeros(round(maxxy(2)-minxy(2))+1,round(maxxy(1)-minxy(1))+1); 
%         im_first_raw_zoom = z; im_second_raw_zoom = z; 
%         im_first_zoom = z; im_second_zoom = z; 
% 
%         im_first_raw_zoom(mask1_small) = im_first_raw(mask1);
%         im_first_raw_zoom = imresize(im_first_raw_zoom,[siz(1),siz(1)]); 
% 
%         im_second_raw_zoom(mask2_small) = im_second_raw(mask2);
%         im_second_raw_zoom = imresize(im_second_raw_zoom,[siz(1),siz(1)]); 
% 
%         im_first_zoom(mask1_small) = im_first(mask1);
%         im_first_zoom = imresize(im_first_zoom,[siz(1),siz(1)]); 
% 
%         im_second_zoom(mask2_small) = im_second(mask2); 
%         im_second_zoom = imresize(im_second_zoom,[siz(1),siz(1)]); 
% 
%         fig = figure; 
%         imshow([im_first_raw, im_first_raw_zoom,...
%                 im_second_raw,im_second_raw_zoom;...
%                 im_first,     im_first_zoom,...
%                 im_second,    im_second_zoom]); 
% 
%         answer = questdlg('Correct points ?', ...
%             'Zoom point', ...
%             'Yes','No','Yes');
%         % Handle response
%         switch answer
%             case 'Yes'
%                 bool_point = false;
%             case 'No'
%                 bool_point = true;
%         end
%         close(fig); 
%     end
%     
%     im_c = zeros(siz(1)*2,siz(2)*2+siz(1)*2,siz(3)); 
%     z = zeros(round(maxxy(2)-minxy(2))+1,round(maxxy(1)-minxy(1))+1); 
%     for ii = 1:siz(3)
%         im_first_raw = cam_first_raw(:,:,ii); 
%         im_second_raw = cam_second_raw(:,:,ii); 
%         im_first = cam_first_satur(:,:,ii); 
%         im_second = cam_second_satur(:,:,ii);
%         
%         %highlight circle
%         im_first_raw(~mask1) = im_first_raw(~mask1)/factor;
%         im_second_raw(~mask2) = im_second_raw(~mask2)/factor;
%         im_first(~mask1) = im_first(~mask1)/factor;
%         im_second(~mask2) = im_second(~mask2)/factor;
%         
%         %zoom
%         im_first_raw_zoom = z; im_second_raw_zoom = z; 
%         im_first_zoom = z; im_second_zoom = z; 
%         
%         im_first_raw_zoom(mask1_small) = im_first_raw(mask1);
%         im_first_raw_zoom = imresize(im_first_raw_zoom,[siz(1),siz(1)]); 
%         
%         im_second_raw_zoom(mask2_small) = im_second_raw(mask2);
%         im_second_raw_zoom = imresize(im_second_raw_zoom,[siz(1),siz(1)]); 
%         
%         im_first_zoom(mask1_small) = im_first(mask1);
%         im_first_zoom = imresize(im_first_zoom,[siz(1),siz(1)]); 
%         
%         im_second_zoom(mask2_small) = im_second(mask2); 
%         im_second_zoom = imresize(im_second_zoom,[siz(1),siz(1)]); 
%         
%         im_c(:,:,ii) = [im_first_raw,im_first_raw_zoom,...
%                         im_second_raw,im_second_raw_zoom;...
%                         im_first,im_first_zoom,...
%                         im_second,im_second_zoom]; 
%     end
%     
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MATCHING STEP - RUN #2 : 
% % From second camera to first camera 
% tic; 
% step1_2_parameters.initial_seed = initial_seed_point_set2.pw;
% % We keep the ROI of the ncorr matching step
% 
% 
% 
% %
% siz = size(cam_second_satur); 
% input1 = cam_first_satur(:,:,1); 
% input2 = zeros(siz(1),siz(2),2,'uint8'); 
% input2(:,:,1) = cam_second_satur(:,:,1); 
% input2(:,:,2) = cam_first_satur(:,:,1); %needed to avoid a bug in ncorr
% h12 = ncorr_dic_rewrited('cam_number',[cam_1,cam_2],...
%     'cam_data_ref',input1,...
%     'cam_data_cur',input2,...
%     'mask',refmask_trial,...
%     'automatic_process',automatic_process,...
%     'step_param',step1_2_parameters,...
%     'base_param',base_parameters);
% 
% % Retreive info for next step 
% ref2_roi = h12.current(1).roi;
% 
% % mapping : format seed point for ncorr2
% U_mapped = h12.data_dic.displacements(1).plot_u_ref_formatted/(subset_spacing+1); 
% V_mapped = h12.data_dic.displacements(1).plot_v_ref_formatted/(subset_spacing+1); 
% 
% initial_seed_point_set2.sw = map_pointcoordinate(initial_seed_point_set1.sw,{U_mapped,V_mapped}); 
% initial_seed_point_set2.pw = map_subset2pixel(initial_seed_point_set2.sw,subset_spacing);
% 
% if isfield(h12,'handles_gui')
%     close(h12.handles_gui.figure); 
% end
% clear('h12');
% 
% elapsedTime = toc; 
% fprintf("--> STEP : Ncorr 1_2 done in %1.1fs\n",elapsedTime); 