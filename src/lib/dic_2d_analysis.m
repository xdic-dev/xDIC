function output = dic_2d_analysis(trial_target, data_path, dic_path, subject_id, phase_id, material, num_pair, nfcond_set, spddxlcond_set, calib_folder_set, ref_trial_id, idx_frame_start, idx_frame_end, frame_jump, showvisu, debug_mode, automatic_process)
    fprintf('Begin of %s\n', mfilename);
    
    % Flatten all combinations (trial_jj, stereopair_kk)
    pairs_trials = combvec(1:num_pair, trial_target)';
    num_combinations = size(pairs_trials, 1);
    
    % Initialize storage for matching results
    matching_results = cell(num_combinations, 1);
    
    % STEP 1: Parallel matching computation
    fprintf('Starting parallel matching computation...\n');
    parfor i = 1:num_combinations
        pair_id = pairs_trials(i, 1);
        trial_id = pairs_trials(i, 2);
        fprintf('Matching: Trial ID: %d, Pair ID: %d\n', trial_id, pair_id);
        
        % Compute matching between stereo pairs
        [refmask_trial_matched, initial_seed_point_set2] = compute_matching(...
            data_path, dic_path, subject_id, material, ...
            sprintf('%03d', trial_id), pair_id, phase_id, ...
            showvisu, idx_frame_start, idx_frame_end, frame_jump, ...
            ref_trial_id, automatic_process);
            
        matching_results{i} = struct(...
            'trial_id', trial_id, ...
            'pair_id', pair_id, ...
            'refmask_trial_matched', refmask_trial_matched, ...
            'initial_seed_point_set2', initial_seed_point_set2);
    end
    
    % STEP 2: Parallel tracking computation
    fprintf('Starting parallel tracking computation...\n');
    % Create combinations for tracking (trial, pair, track_id)
    tracking_combinations = [];
    for i = 1:num_combinations
        tracking_combinations = [tracking_combinations; ...
            pairs_trials(i,:), 1; ... % track_id 1 for first camera
            pairs_trials(i,:), 2];    % track_id 2 for second camera
    end
    
    num_tracking_tasks = size(tracking_combinations, 1);
    tracking_results = cell(num_tracking_tasks, 1);
    
    parfor i = 1:num_tracking_tasks
        trial_id = tracking_combinations(i, 2);
        pair_id = tracking_combinations(i, 1);
        track_id = tracking_combinations(i, 3);
        
        fprintf('Tracking: Trial ID: %d, Pair ID: %d, Track ID: %d\n', ...
            trial_id, pair_id, track_id);
        
        % Find matching results for this trial/pair
        matching_idx = -1;
        for j = 1:numel(matching_results)
            if matching_results{j}.trial_id == trial_id && matching_results{j}.pair_id == pair_id
                matching_idx = j;
                break;
            end
        end
        match_data = matching_results{matching_idx};
        
        % Get appropriate mask and seed point based on track_id
        if track_id == 1
            refmask = match_data.refmask_trial_matched;
            initial_seed_point = match_data.initial_seed_point_set2;
        else
            refmask = match_data.refmask_trial_matched;
            initial_seed_point = match_data.initial_seed_point_set2;
        end
        
        % Compute tracking
        success = compute_tracking(...
            data_path, dic_path, subject_id, material, ...
            sprintf('%03d', trial_id), pair_id, phase_id, ...
            showvisu, idx_frame_start, idx_frame_end, frame_jump, ...
            ref_trial_id, automatic_process, track_id, ...
            refmask, initial_seed_point);
            
        tracking_results{i} = struct(...
            'trial_id', trial_id, ...
            'pair_id', pair_id, ...
            'track_id', track_id, ...
            'success', success);
    end
    
    % Combine results
    output = cell(num_combinations, 1);
    for i = 1:num_combinations
        trial_id = pairs_trials(i, 2);
        pair_id = pairs_trials(i, 1);
        
        % Get tracking results for this trial/pair
        track_success = true;
        for j = 1:numel(tracking_results)
            if tracking_results{j}.trial_id == trial_id && tracking_results{j}.pair_id == pair_id
                track_success = track_success && tracking_results{j}.success;
            end
        end
        
        % Format output as before
        output{i} = {sprintf('%03d', trial_id), ...
                    num2str(pair_id), ...
                    num2str(track_success)};
    end
    
    fprintf('End of %s\n', mfilename);
end
