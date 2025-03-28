function file_reconstruction = dic_3d_reconstruction(output_2d, trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process)
    fprintf('Begin of %s\n', mfilename);

    % calib folder path
    calibPath = fullfile(dic_path, subject_id, "calib", calib_folder_set);

    % parfor over the trials
    file_reconstruction = cell(size(trial_target, 1), 1);

    parfor i = 1:size(trial_target, 1)
        trial_id = trial_target(i);
        fprintf('Trial ID: %d\n', trial_id);

        file_reconstruction{i} = stepE_3DReconstruction('basePath', output_2d{i}{1}, ...
            'calibPath', calibPath, ...
            'showvisu', showvisu, ...
            'pairOrder', output_2d{i}{2}, ...
            'pairForced', output_2d{i}{3}, ...
            'title_ax', sprintf("%s-%s-%s-%s", subject_id, material, sprintf("%03d", trial_id), phase_id));
    end
    fprintf('End of %s\n', mfilename);
end