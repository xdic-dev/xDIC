function perform_tracking(config, base_parameters, step_parameters)
% PERFORM_TRACKING Helper function to perform tracking for a single camera
%   This function executes the tracking analysis for a single camera using the
%   provided configuration and parameters.
%
% Parameters:
%   config: Structure containing tracking configuration for one camera
%       - cam_number: Camera number identifier
%       - ref_data: Reference image data
%       - cur_data: Current image data for tracking
%       - mask: Mask for the region of interest
%       - initial_seed: Initial seed point for tracking
%       - step_number: Step number identifier
%   base_parameters: Structure containing base parameters
%   step_parameters: Structure containing step-specific parameters

tic;

% Set initial seed for this tracking step
step_parameters.initial_seed = config{1}.initial_seed;

% Perform tracking analysis
[h, file_logic] = ncorr_dic_rewrited(...
    'cam_number', config{1}.cam_number,...
    'cam_data_ref', config{1}.ref_data,...
    'cam_data_cur', config{1}.cur_data,...
    'mask', config{1}.mask,...
    'automatic_process', base_parameters.automatic_process,...
    'step_param', step_parameters,...
    'base_param', base_parameters);

% Clean up
if ~file_logic
    close(h.handles_gui.figure);
end
clear('h');

% Report timing
elapsedTime = toc;
fprintf('--> STEP: Ncorr %d done in %1.1fs\n', config{1}.step_number, elapsedTime);
end
