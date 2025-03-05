
% Frictional conditions
frictional_conditions = ["glass", "coating", "coating_oil"];

% Global parameters
num_pair = 2; % Number of stereopairs
robot_sample_freq = 1e3; % robot sampling frequency (Hz)
vid_sample_freq = 50; % image sampling frequency (Hz)

% Path definitions
base_path = "./"; % base path
data_path = base_path + "data"; % location of the video files and protocol files
DIC_path = base_path + "analysis"; % location of the output data from DIC

% processing flags
automatic_process = true; % automatic processing flag
