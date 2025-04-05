function [initial_seed_point_set1, refmask_trial] = compute_initial_seed_points(cam_first_satur, cam_second_satur)
    % This function computes initial seed points for DIC analysis
    % It is extracted from stepD_2DDIC.m to support parallel processing
    
    % Initialize output
    initial_seed_point_set1 = struct();
    
    % Get size of images
    siz = size(cam_first_satur);
    
    % Create initial mask (full image)
    refmask_trial = true(siz(1:2));
    
    % Set initial seed point in the center of the image
    % This is a simple approach - you might want to make this more sophisticated
    initial_seed_point_set1.sw = [siz(2)/2, siz(1)/2];
    initial_seed_point_set1.pw = initial_seed_point_set1.sw;
end
