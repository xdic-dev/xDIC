
function varargout = ncorr_dic_rewrited(varargin)

disp('Ncorr analysis launch');
%% Parse input arguments
p = inputParser;
p.addParameter('cam_number',[]);
p.addParameter('cam_data_ref',[]);
p.addParameter('cam_data_cur',[]);
p.addParameter('mask',[]);
p.addParameter('automatic_process',[]);
p.addParameter('step_param',[]);
p.addParameter('base_param',[]);

p.parse(varargin{:});
cam_number = p.Results.cam_number;
cam_data_ref = p.Results.cam_data_ref;
cam_data_cur = p.Results.cam_data_cur;
mask = p.Results.mask;
automatic_process = p.Results.automatic_process;
step_param = p.Results.step_param;
base_param = p.Results.base_param;

%% 
filename = fullfile(base_param.outputPath, "ncorr"+sprintf("%d",cam_number)+".mat");

%% DIC process
if ~exist(filename,'file') %make sure the file doesn't already exist before doing the DIC process
file_logic = false; 
if (isfield(step_param, "initial_seed"))
    initial_seed = step_param.initial_seed;
else
     initial_seed = [];
end
showGui = ~automatic_process; 
handles = ncorr_abr(showGui,initial_seed);
% handles = ncorr;
handles.set_ref(cam_data_ref);
handles.set_cur(Myarray2cell(cam_data_cur));

handles.set_roi_ref(mask);


handles.set_dic_parameters_manually( ...
    step_param.type, ...
    step_param.radius, ...
    step_param.spacing, ...
    step_param.cutoff_diffnorm, ...
    step_param.cutoff_iteration,...
    step_param.total_threads, ...
    step_param.stepanalysis_params)

% refresh the handle
handles.abr_refresh()

% trigger DIC
handles.callback_topmenu_dic();

% trigger Displacement
handles.callback_topmenu_formatdisp();


%% save ncorr results
save_dicdata(handles, filename); 
% close(handles.handles_gui.figure);    
else 
file_logic = true; 
handles = load(filename);
handles.current = handles.current_save; handles = rmfield(handles,'current_save');
handles.data_dic = handles.data_dic_save; handles = rmfield(handles,'data_dic_save');
handles.reference = handles.reference_save; handles = rmfield(handles,'reference_save');
end
%% output argument
if nargout >= 1
    varargout{1} = handles; 
end
if nargout == 2
    varargout{2} = file_logic; 
end

end














