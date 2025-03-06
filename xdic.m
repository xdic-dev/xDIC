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

% Add path
addpath(genpath(pwd));

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


