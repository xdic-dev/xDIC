function varargout = import_raw_vid(basePath,varargin)
% IMPORT_RAW_VID imports video pairs files from MP4 extension to the MATLAB workspace.
% The function takes the following input arguments:
%   import_raw_vid(basePath, 'subject', subject_value, 'trial', trial_value, 'stereopair', stereopair_value)
%
% Input:
%   - basePath: The base path where the video files are located.
%   - subject: The subject identifier.
%   - trial: The trial identifier.
%   - stereopair: The stereo pair identifier.
%
% Output:
%   - varargout{1}: The first video file imported as a MATLAB variable.
%   - varargout{2}: The second video file imported as a MATLAB variable.

% Parse input arguments
p = inputParser;
p.addRequired('basePath', @isText);
p.addParameter('subject',[], @isText);
p.addParameter('material',[], @isText);
p.addParameter('trial',[], @isText);
p.addParameter('stereopair',[], @isnumeric);
p.addParameter('frameStart',[], @isnumeric);
p.addParameter('frameEnd',[], @isnumeric);
p.addParameter('framejump',[], @isnumeric);

p.parse(basePath, varargin{:});
subject = p.Results.subject;
trial = p.Results.trial;
material = p.Results.material;
stereopair = p.Results.stereopair;
frameStart = p.Results.frameStart;
frameEnd = p.Results.frameEnd;
framejump = p.Results.framejump;

trialname = trial; %subject+"_"+material+"_speckles_"+trial;

% Identify cameras based on stereopair value
switch stereopair
    case 1
        % Pair on the radial side of the index fingertip:
        % The first camera is the one most deported from the front.
        cam_first = 1;
        cam_second = 2;
    case 2
        % Pair on the ulnar side of the index fingertip.
        cam_first = 4;
        cam_second = 3;
    otherwise
        error('Error: Unrecognized stereopair value.');
end

% Find and import the video files
vid_file_pattern_first = fullfile(basePath,subject,"speckles",material,"vid",sprintf('*_%s*_cam_%d*.mp4', trialname, cam_first));
S = dir(vid_file_pattern_first);
if isempty(S)
    sprintf('*_%s*_cam_%d*.mp4', trialname, cam_first)
    error('Error: First video file not found.');
end
vid_first = readvid(fullfile(S.folder, S.name),(frameEnd-frameStart+1)/framejump,frameStart-1,framejump);

vid_file_pattern_second = fullfile(basePath,subject,"speckles",material,"vid",sprintf('*_%s*_cam_%d*.mp4', trialname, cam_second));
S = dir(vid_file_pattern_second);
if isempty(S)
    error('Error: Second video file not found.');
end
vid_second = readvid(fullfile(S.folder, S.name),(frameEnd-frameStart+1)/framejump,frameStart-1,framejump);

% Assign output variables
varargout{1} = vid_first;
varargout{2} = vid_second;
if nargout >= 2
    varargout{3} = cam_first; 
    varargout{4} = cam_second; 
end
end

function output = isText(myVar)
% isText checks if the input variable is a text type.
% The function returns true if the variable is of string, character array, or cell array of strings type; otherwise, it returns false.
%
% Input:
%   - myVar: The variable to be checked.
%
% Output:
%   - output: A logical value indicating if the variable is a text type.

    output = isstring(myVar) || ischar(myVar) || iscellstr(myVar);
end