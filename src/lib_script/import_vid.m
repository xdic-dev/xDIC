function varargout = import_vid(basePath,varargin)
% IMPORT_VID imports video files from the specified basePath according to the given parameters.
%
% Input:
%   - basePath: The base path where the video files are located.
%   - subject: The subject identifier.
%   - trial: The trial identifier.
%   - stereopair: The stereo pair identifier.
%   - phase: The type of phase for frame selection.
%   - framejump: The frame jump value for subsampling.
%
% Output:
%   - varargout{1}: The frames from the first camera as a MATLAB variable.
%   - varargout{2}: The frames from the second camera as a MATLAB variable.

% Parse input arguments
p = inputParser;
p.addRequired('basePath', @isText);
p.addParameter('subject',[],@isText);
p.addParameter('material',[],@isText);
p.addParameter('trial',[],@isText);
p.addParameter('stereopair',[],@isnumeric);
p.addParameter('idxstart_set',[],@isnumeric);
p.addParameter('idxend_set',[],@isnumeric);
p.addParameter('phase',[],@isText);
p.addParameter('framejump',1,@isnumeric);

p.parse(basePath, varargin{:});
subject = p.Results.subject;
material = p.Results.material;
trial = p.Results.trial;
stereopair = p.Results.stereopair;
phase = p.Results.phase;
idxstart_set = p.Results.idxstart_set;
idxend_set = p.Results.idxend_set;
framejump = p.Results.framejump;

% Validate input arguments
validateInputArguments(basePath, subject, trial, stereopair, phase, framejump);

% Load protocol file
p = fullfile(basePath,'rawdata',subject,'speckles',material,'protocol','*.mat');
S = dir(p);
protocolPath = fullfile(S.folder,S.name);
validateProtocolFile(protocolPath);
protocol = load(protocolPath);

% Select frames according to protocol
protocol = protocol.cond;
trialNum = str2double(trial); 
dst = protocol.table{trialNum,5}; % distance in mm
spd = protocol.table{trialNum,4}; % velocity of phase in mm/s
fps = 50; % frame per second in Hz 

NbrFrLoading = protocol.dur(3)/1e3*fps; 
NbrFrSlide = dst/spd*fps; 
NbrFrRelax = protocol.dur(5)/1e3*fps; 

% Calculate frame ranges based on the specified phase
[frameRangeStart, frameRangeEnd] = calculateFrameRanges(phase, NbrFrLoading, NbrFrSlide, NbrFrRelax);
if ~isempty(idxstart_set)
    frameRangeStart = idxstart_set; 
    frameRangeEnd = idxend_set;
end
% Import raw video files
[raw_cam_first, raw_cam_second, cam1, cam2] = import_raw_vid(fullfile(basePath,'rawdata'),...
                                                   'subject',subject,...
                                                   'material',material,...
                                                   'trial',trial,...
                                                   'stereopair',stereopair,...
                                                   'frameStart',frameRangeStart,...
                                                   'frameEnd',frameRangeEnd,...
                                                   'framejump',framejump...
                                               );
% Assign output variables
varargout{1} = raw_cam_first;
varargout{2} = raw_cam_second;
if nargout >= 2 
    varargout{3} = cam1;  
    varargout{4} = cam2; 
end
end

function validateProtocolFile(protocolPath)
% Validate the protocol file and throw an error if it fails to load or is missing

if ~isfile(protocolPath)
    error('Protocol file is missing. Please make sure the protocol file exists in the specified path.');
end

try
    protocol = load(protocolPath);
catch
    error('Failed to load the protocol file. Please make sure the protocol file is valid.');
end

% Check if the required fields exist in the protocol structure
if ~isfield(protocol, 'cond') || ~isfield(protocol.cond, 'table') || ~isfield(protocol.cond, 'dur')
    error('Invalid protocol file. Please make sure the protocol file contains the required fields.');
end
end

function validateInputArguments(basePath, subject, trial, stereopair, phase, framejump)
% Validate input arguments and throw errors if any requirements are not met

% Check if basePath is provided and exists
if isempty(basePath) || ~isfolder(basePath)
    error('Invalid or missing basePath. Please provide a valid base path.');
end

% Check if subject, trial, stereopair, and phase are provided
if isempty(subject) || isempty(trial) || isempty(stereopair) || isempty(phase)
    error('Missing required input arguments. Please provide subject, trial, stereopair, and phase.');
end

% Check if framejump is provided and a positive integer
if isempty(framejump) || ~isnumeric(framejump) || framejump <= 0 || mod(framejump, 1) ~= 0
    error('Invalid framejump value. Please provide a positive integer value for framejump.');
end
end

