
% Frictional conditions
frictional_conditions = ["glass", "coating", "coating_oil"];

% Global parameters
num_pair = 2; % Number of stereopairs
robot_sample_freq = 1e3; % robot sampling frequency (Hz)
vid_sample_freq = 50; % image sampling frequency (Hz)

% Path definitions
base_path = "E:\MATLAB\multiDIC"; % base path
data_path = base_path + "\example_data"; % location of the video files and protocol files
dic_path = base_path + "\xDIC\analysis"; % location of the output data from DIC

% processing flags
automatic_process = true; % automatic processing flag
parallel_processing = true; % parallel processing flag
