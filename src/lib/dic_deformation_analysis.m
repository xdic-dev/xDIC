function dic_deformation_analysis(file_reconstruction, output_2d, trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process)
    fprintf('Begin of %s\n', mfilename);

    % parfor over the trials
    parfor i = 1:size(trial_target, 1)
        trial_id = trial_target(i);
        fprintf('Trial ID: %d\n', trial_id);
        tic;

        % Deformation analysis
        stepF_Deformation('basePath', output_2d{i}{1}, ...
            'target_file', file_reconstruction{i});
            
        elapsed_time = toc; 
        fprintf("DIC Deformation Analysis done in %1.1fs\n", elapsed_time);
    end
        
