function dic_analysis(data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process)
    %trial_target = search_trial2target(data_path, subject_id, phase_id, material, nfcond_set, spddxlcond_set);
    trial_target = [5, 6];
    fprintf('Trial target set: [%s]\n', num2str(trial_target));

    % STEP D: 2D-DIC
    % 2D-DIC analysis
    tic;
    output = dic_2d_analysis(trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);
    elapsed_time = toc;
    fprintf("Dic 2D Analysis done in %lf s\n", elapsed_time);

    % output_2d
    output_2d = cell(size(output, 1)/num_pair, 1);
    j = 1;
    for i = 1:size(output, 1)
        if mod(i, num_pair) == 1
            output_2d{j} = {output{i}(1), [str2num(output{i}(2)), str2num(output{i}(3))], str2num(output{i}(4))};
            j = j + 1;
        end
    end

    % STEP E: 3D Reconstruction
    % 3D reconstruction
    tic;
    file_reconstruction = dic_3d_reconstruction(output_2d, trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);
    elapsed_time = toc; 
    disp(file_reconstruction);
    fprintf("DIC 3D Reconstruction done in %1.1fs\n", elapsed_time);

    % STEP F: Deformation analysis
    % Deformation analysis
    tic;
    dic_deformation_analysis(file_reconstruction, output_2d, trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process);
    elapsed_time = toc; 
    fprintf("DIC Deformation Analysis done in %1.1fs\n", elapsed_time);
 
end
