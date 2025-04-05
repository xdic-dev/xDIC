function [matching_results] = dic_matching(prep_params)
% Get required parameters
base_parameters = prep_params.base_parameters;
step1_2_parameters = prep_params.step1_2_parameters;
cam_data = prep_params.cam_data;

% Define file paths
roifile = fullfile(base_parameters.baseResultPath, base_parameters.subject, ...
    base_parameters.material, sprintf("REF_MASK_%s_%s_pair%d.mat", ...
    prep_params.reftrial, base_parameters.phase, base_parameters.stereopair));

matchingfile = fullfile(prep_params.outputPath, ...
    sprintf("MATCHING2%s_pair%d.mat", prep_params.reftrial, base_parameters.stereopair));

seedfile = fullfile(base_parameters.baseResultPath, base_parameters.subject, ...
    base_parameters.material, sprintf("REF_SEED_%s_%s_pair%d.mat", ...
    prep_params.reftrial, base_parameters.phase, base_parameters.stereopair));

%% ROI and Seed initialization
% Draw ROI if needed
if ~exist(roifile, 'file')
    draw_ref_roi(base_parameters);
end

% Perform matching if needed
if ~exist(matchingfile, 'file')
    ncorr_matching2ref(cam_data.first_satur(:,:,1), base_parameters, step1_2_parameters);
end

% Load matching results
matching = load(matchingfile);
refmask_REF = matching.reference_save(1).roi.mask;
refmask_trial = matching.current_save(1).roi.mask;
fprintf('--> STEP: ROI loaded and formatted\n');

% Handle seed points
if ~exist(seedfile, 'file')
    draw_ref_seed(base_parameters, 'roi', refmask_REF);
end

% Load and map seed points
seed_point = load(seedfile);
ref_seed_point.pw = seed_point.seed_point;
ref_seed_point.sw = map_pixel2subset(ref_seed_point.pw, step1_2_parameters.spacing);

% Map coordinates
U_mapped = matching.data_dic_save.displacements(1).plot_u_ref_formatted/(step1_2_parameters.spacing + 1);
V_mapped = matching.data_dic_save.displacements(1).plot_v_ref_formatted/(step1_2_parameters.spacing + 1);

initial_seed_point_set1.sw = map_pointcoordinate(ref_seed_point.sw, {U_mapped, V_mapped});
initial_seed_point_set1.pw = map_subset2pixel(initial_seed_point_set1.sw, step1_2_parameters.spacing);

fprintf('--> STEP: SEED loaded and formatted\n');

%% Perform matching analysis
tic;
step1_2_parameters.initial_seed = initial_seed_point_set1.pw;

siz = size(cam_data.first);
input1 = cam_data.first_satur(:,:,1);
input2 = zeros(siz(1), siz(2), 2, 'uint8');
input2(:,:,1) = cam_data.second_satur(:,:,1);
input2(:,:,2) = cam_data.first_satur(:,:,1);

[h12, file_logic] = ncorr_dic_rewrited('cam_number', [cam_data.cam_1, cam_data.cam_2],...
    'cam_data_ref', input1,...
    'cam_data_cur', input2,...
    'mask', refmask_trial,...
    'automatic_process', base_parameters.automatic_process,...
    'step_param', step1_2_parameters,...
    'base_param', base_parameters);

% Get results for next step
refmask_trial_matched = h12.current(1).roi.mask;
U_mapped = h12.data_dic.displacements(1).plot_u_ref_formatted/(step1_2_parameters.spacing + 1);
V_mapped = h12.data_dic.displacements(1).plot_v_ref_formatted/(step1_2_parameters.spacing + 1);

initial_seed_point_set2.sw = map_pointcoordinate(initial_seed_point_set1.sw, {U_mapped, V_mapped});
initial_seed_point_set2.pw = map_subset2pixel(initial_seed_point_set2.sw, step1_2_parameters.spacing);

if ~file_logic
    close(h12.handles_gui.figure);
end
clear('h12');

elapsedTime = toc;
fprintf("--> STEP: Ncorr matching 1-2 done in %1.1fs\n", elapsedTime);

%% Prepare output structure
matching_results = struct(...
    'refmask_trial', refmask_trial,...
    'refmask_trial_matched', refmask_trial_matched,...
    'initial_seed_point_set1', initial_seed_point_set1,...
    'initial_seed_point_set2', initial_seed_point_set2);
end
