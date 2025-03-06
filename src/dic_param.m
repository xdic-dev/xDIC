% STEP D, E, F: 2D-DIC, 3D Reconstruction, Deformation analysis
subject_id = "S09"; % Subject identifier
phase_id = "loading";  % Phase identifier
material_id = 2; % Material identifier
nfcond_set = [5]; % Number of conditions
spddxlcond_set = [0.04, 0.08]; % Speed conditions

% Calibration settings
calib_folder_set = "2";

% Reference trial number
ref_trial_id = 5;

% Frame settings
idx_frame_start = 1;
idx_frame_end = 150;
frame_jump = 1;

% Visualization settings
showvisu = 0; % Boolean for visualization
debug_mode = 0; % Debug mode flag
