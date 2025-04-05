
function varargout = stepD_2DDIC(varargin)
% STEPD_2DDIC Main function for Digital Image Correlation analysis
%   This function orchestrates the DIC analysis pipeline by:
%   1. Preparing data and parameters (using dic_preparation)
%   2. Performing matching between image pairs (using dic_matching)
%   3. Executing tracking analysis for both cameras
%
% Parameters:
%   varargin: Variable input arguments including:
%     - baseDataPath: Base path for data files
%     - baseResultPath: Base path for results
%     - subject: Subject identifier
%     - material: Material type
%     - trial: Trial number
%     - stereopair: Stereo pair number
%     - phase: Analysis phase
%     - jump: Frame jump size
%     - idxstart_set: Start index
%     - idxend_set: End index
%     - showvisu: Visualization flag
%     - savedata: Save data flag
%     - reftrial_setmanual: Manual reference trial
%     - param_filt_im: Image filtering parameters
%     - automatic_process: Automatic processing flag
%
% Returns:
%   varargout: Optional output structure containing:
%     - prep_params: Parameters and data from preparation step
%     - matching_results: Results from the matching analysis
%
% See also:
%   dic_preparation, dic_matching, perform_tracking

% Print analysis start message
fprintf('-------------------------------------------\n');
fprintf('-------------------------------------------\n');
fprintf('Digital Image Correlation analysis launch\n');

% Prepare data and parameters
prep_params = dic_preparation(varargin{:});

% Perform matching analysis
matching_results = dic_matching(prep_params);

% Execute tracking analysis for both cameras
configs = prep_params.tracking_configs;
configs{1}.mask = matching_results.refmask_trial;
configs{1}.initial_seed = matching_results.initial_seed_point_set1.pw;
configs{2}.mask = matching_results.refmask_trial_matched;
configs{2}.initial_seed = matching_results.initial_seed_point_set2.pw;

for config = configs
    perform_tracking(config, prep_params.base_parameters, prep_params.step1_parameters);
end

% Return any output arguments if needed
if nargout > 0
    varargout{1} = struct(...
        'prep_params', prep_params,...
        'matching_results', matching_results);
end

% Print completion message
fprintf('-------------------------------------------\n');
fprintf('Digital Image Correlation analysis complete\n');
fprintf('-------------------------------------------\n');
end