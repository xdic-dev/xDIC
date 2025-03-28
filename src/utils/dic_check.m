% parallel checking
if parallel_processing
    disp('Parallel processing is enabled');
    %delete(gcp('nocreate'));
    
    disp('Checking ROI References...');
    for pair_i = 1:num_pair
        roifile = fullfile(dic_path, subject_id, material, "REF_MASK_"+sprintf("%03d", ref_trial_id)+"_"+phase_id+"_pair"+num2str(pair_i)+".mat");
        if ~exist(roifile,'file')    
            fprintf("This file %s must exist to be able to run DIC analysis in parallel.\n", roifile)
            fprintf("This file come from Drawing Reference ROI process you can run it using this command or just copy from previous experience into the right place.\n")
            error('Error: ROI file for the stereopair %d from the Reference Trial %d for the subject %s not found.\n', pair_i, ref_trial_id, subject_id);
        end
    end

    disp('Checking SEED References...');
    for pair_i = 1:num_pair
        seedfile = fullfile(dic_path, subject_id, material, "REF_SEED_"+sprintf("%03d", ref_trial_id)+"_"+phase_id+"_pair"+num2str(pair_i)+".mat");
        if ~exist(seedfile,'file')    
            fprintf("This file %s must exist to be able to run DIC analysis in parallel.\n", seedfile)
            fprintf("This file come from Selecting the seed of ROI process you can run it using this command or just copy from previous experience into the right place.\n")
            error('Error: SEED file for the stereopair %d from the Reference Trial %d for the subject %s not found.\n', pair_i, ref_trial_id, subject_id);
        end
    end

end

protocol_path = fullfile(fullfile(data_path,"rawdata",subject_id,"speckles", material,"protocol",sprintf('*.mat')));
if isempty(dir(protocol_path))
    fprintf("These files %s must exist to be able to run DIC analysis.\n", protocol_path)
    fprintf("These files come from xxx process you can run it using this command or just copy from previous experience into the right place.\n")
    error('Error: Protocol not found.\n');
end

calib_path = fullfile(fullfile(dic_path,subject_id,"calib", calib_folder_set,sprintf('*cam_*.mat')));
if isempty(dir(calib_path))
    fprintf("These files %s must exist to be able to run DIC analysis useful for 3D reconstruction.\n", calib_path)
    fprintf("These files come from xxx process you can run it using this command or just copy from previous experience into the right place.\n")
    error('Error: Calib files not found.\n');
end

disp('...done Checking.');