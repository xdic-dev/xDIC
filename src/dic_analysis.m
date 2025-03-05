function dic_analysis(data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process)
    %trial_target = search_trial2target(baseDataPath, subject, phase, material, nfcond_set, spddxlcond_set);
    trial_target = [5, 7, 9, 10];
    fprintf('Trial target set: [%s]\n', num2str(trial_target));

    % STEP D: 2D-DIC
    % 2D-DIC analysis
    dic_2d_analysis(trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);

    % STEP E: 3D Reconstruction
    % 3D reconstruction
    %dic_3d_reconstruction(trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);

    % STEP F: Deformation analysis
    % Deformation analysis
    %dic_deformation_analysis(trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);
 
end