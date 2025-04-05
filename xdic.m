        %%% Main script for FINGERTIP 3D RECONSTRUCTION using DIC
% This script is the main script for the 3D reconstruction of the fingertip
% using DIC. It calls the analysis function and sets the parameters for the
% analysis.

% Run the following command in the terminal to execute this script:
% matlab -nodisplay -nojvm -nosplash -nodesktop -r "run('xdic.m'); exit;"

clear all
close all
clc

%
disp('FINGERTIP 3D RECONSTRUCTION using DIC');
disp('--------------------------------------');
disp('DIC analysis for the fingertip');

% Add paths
addpath(genpath(pwd));

% Add ncorr path and cd to it for initialization
ncorr_path = fullfile(pwd, 'toolbox', 'MultiDIC-master', 'lib_ext', 'ncorr_2D_matlab-master');
addpath(ncorr_path);
current_dir = pwd;
cd(ncorr_path);

% Ensure workers have access to the code
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool('local');

% Attach all necessary files to the parallel pool
lib_files = {...
    'compute_initial_seed_points.m', ...
    'compute_matching.m', ...
    'compute_tracking.m', ...
    'dic_2d_analysis.m', ...
    'dic_3d_reconstruction.m', ...
    'dic_deformation_analysis.m', ...
    'stepD_2DDIC.m', ...
    'stepE_3DReconstruction.m', ...
    'stepF_Deformation.m' ...
};

for i = 1:length(lib_files)
    lib_files{i} = fullfile(current_dir, 'src', 'lib', lib_files{i});
end

% Add ncorr files
ncorr_files = dir(fullfile(ncorr_path, '*.m'));
for i = 1:length(ncorr_files)
    lib_files{end+1} = fullfile(ncorr_files(i).folder, ncorr_files(i).name);
end

% Attach files and update path on workers
addAttachedFiles(gcp, lib_files);
parfevalOnAll(gcp, @(p) cd(p), 1, ncorr_path);
parfevalOnAll(gcp, @(p) addpath(p, genpath(p)), 1, ncorr_path);

% Return to original directory
cd(current_dir);

% Load configurations
global_param;
dic_param;
update_variables;

% print parameters
disp('Parameters:');
disp('-----------');
fprintf('Subject: %s, Phase: %s, Material: %s, Stereo Pairs: %d\n', subject_id, phase_id, material, num_pair);
fprintf('Reference trial number: %d\n', ref_trial_id);
fprintf('Frame: %d to %d, jump= %d\n', idx_frame_start, idx_frame_end, frame_jump);
fprintf('Show visualization: %d, Debug mode: %d, Automatic process: %d\n', showvisu, debug_mode, automatic_process);
fprintf('\n');

% checking
disp('Checking the data and protocol...');
dic_check;

% call analysis function
dic_analysis(data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);

% end
disp('End of script');


