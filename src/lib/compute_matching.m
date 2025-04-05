function [refmask_trial_matched, initial_seed_point_set2] = compute_matching(data_path, dic_path, subject_id, material, trial, stereopair, phase, showvisu, idx_frame_start, idx_frame_end, frame_jump, ref_trial_id, automatic_process)
    % This function handles the matching step between stereo pairs
    % It is extracted from stepD_2DDIC.m to enable parallel processing
    
    % Get parameters from original function
    TRUE_FPS = 50;
    if str2double(cell2mat(regexp(subject_id, '\d+', 'match'))) < 8
        LIMIT_GRAYSCALE = 70;
    else
        LIMIT_GRAYSCALE = 100;
    end
    
    % Ncorr DIC parameters for matching
    analysis_direction = 'regular';
    subset_radius_ncorr_matching = 60;
    subset_spacing = 10;
    cutoff_matching = 1e-5;
    number_iteration_solver = 100;
    number_threads = 1;
    high_strain_analysis = 1;
    seed_propagation = 'seed';
    auto_ref_change = 1;
    step_ref_change = 10;

    % Set up parameters struct
    step1_2_parameters = struct(...
        'type', analysis_direction,...
        'radius', subset_radius_ncorr_matching,...
        'spacing', subset_spacing,...
        'cutoff_diffnorm', cutoff_matching,...
        'num_threads', number_threads,...
        'enabled_stepanalysis', high_strain_analysis,...
        'seedpropagation', seed_propagation,...
        'autorefchange', auto_ref_change,...
        'steprefchange', step_ref_change);

    base_parameters = struct(...
        'num_iterations', number_iteration_solver);

    % Read video data and process images
    [cam_first_raw, cam_second_raw, cam_1, cam_2] = import_vid(data_path,...
        'subject', subject_id,...
        'material', material,...
        'trial', trial,...
        'stereopair', stereopair,...
        'phase', phase,...
        'idxstart_set', idx_frame_start,...
        'idxend_set', idx_frame_end,...
        'framejump', frame_jump);

    if strcmp(phase, "slide1")
        cam_first_raw = cam_first_raw(:,:,1:end/2+5);
        cam_second_raw = cam_second_raw(:,:,1:end/2+5);
    end

    % Image saturation
    cam_first_satur = satur(cam_first_raw, 'level', LIMIT_GRAYSCALE);
    cam_second_satur = satur(cam_second_raw, 'level', LIMIT_GRAYSCALE);

    % Get initial seed points (this would need to be passed or computed)
    [initial_seed_point_set1, refmask_trial] = compute_initial_seed_points(cam_first_satur, cam_second_satur);
    
    % Create output directory
    outputPath = fullfile(dic_path, subject_id, material, trial, phase);
    if ~exist(outputPath, 'dir')
        mkdir(outputPath);
    end
    base_parameters.outputPath = outputPath;

    % MATCHING STEP
    step1_2_parameters.initial_seed = initial_seed_point_set1.pw;
    siz = size(cam_first_raw);
    input1 = cam_first_satur(:,:,1);
    input2 = zeros(siz(1), siz(2), 2, 'uint8');
    input2(:,:,1) = cam_second_satur(:,:,1);
    input2(:,:,2) = cam_first_satur(:,:,1);

    % Save info before ncorr analysis
    savefileName = fullfile(outputPath, sprintf('dic_info_data_target_pair%d.mat', stereopair));
    actual_fps_meas = TRUE_FPS/frame_jump;
    idxframe = idx_frame_start:frame_jump:idx_frame_end;
    save(savefileName, 'actual_fps_meas', 'idxframe');

    % Run ncorr analysis
    [h12, file_logic] = ncorr_dic_rewrited('cam_number', [cam_1, cam_2],...
        'cam_data_ref', input1,...
        'cam_data_cur', input2,...
        'mask', refmask_trial,...
        'automatic_process', automatic_process,...
        'step_param', step1_2_parameters,...
        'base_param', base_parameters);

    % Get results
    refmask_trial_matched = h12.current(1).roi.mask;
    U_mapped = h12.data_dic.displacements(1).plot_u_ref_formatted/(subset_spacing+1);
    V_mapped = h12.data_dic.displacements(1).plot_v_ref_formatted/(subset_spacing+1);
    
    initial_seed_point_set2.sw = map_pointcoordinate(initial_seed_point_set1.sw, {U_mapped, V_mapped});
    initial_seed_point_set2.pw = map_subset2pixel(initial_seed_point_set2.sw, subset_spacing);

    if ~file_logic
        close(h12.handles_gui.figure);
    end
    clear('h12');
end
