function signals = smp_import(fpath,varargin)
%SMP_IMPORT Extract and filter SMP force, position and acceleration data.
% Adapted from at_import from libactivetouch
% Data are assumed to be acquired at 1000 Hz.
%
% Syntax: signals = smp_import(fpath,runfilter)
%
% Inputs:
%   fpath      path to the SMP data file
%   runfilter  [opt] true = filter data (default); false = use raw data
%
% Outputs:
%   signals   filtered force, acc and COP data, dim = (nbsamples,nsigs)
%             rows: one row for each time step (data acquired at 200 Hz)
%             columns: time,frames,gf,lf,lfv,lfh  (time starts at 0, frames start at 1)
%                      fx,fy,fz,ft,tx,ty,tz,      for CAM side ATI [N,Nm]
%                      fx,fy,fz,ft,tx,ty,tz,      for LED side ATI [N,Nm]
%                      copx,copy,                 for CAM side ATI [mm]
%                      copx,copy,                 for LED side ATI [mm]
%                      ax,ay,az                   for manipulandum [m/s^2]
%
% GF = mean normal force on sensors (along Z)
% LFv = sum of vertical tangential force on sensors (along Y)
% LFh = sum of horizontal tangential force on sensors (along X)
% LF = sqrt(LFv^2 + LFh^2)

% Inputs
if(nargin < 2), varargin{1} = false; varargin{2} = 100; end
runfilter = varargin{1};
freqFilt = varargin{2};
% Constants
freqAcq = 1000;      % acquisition frequency [Hz] 
imAcq = 100;         % images acquisition rate [fps]

% Tools
[B,A] = butter(4,freqFilt/(freqAcq/2)); % low-pass filter for force sigs

% Load data
data = importdata(fpath);
data = data.data; %comment for "January to June 2022" data from SMP

% Extract data
cmd_robot = data(:,1:4);  
position_robot_raw = data(:,5:7);
angle_robot = data(:,4);
forces_mini_raw = data(:,9:14);
forces_nano_raw = data(:,15:20);
moisture = data(:,21);
acc = data(:,22:24);
time = ((0:(size(data,1)-1))/freqAcq)';

% Filter force signals
if(runfilter)
    forces_mini = filtfilt(B,A,forces_mini_raw);
    forces_nano = filtfilt(B,A,forces_nano_raw);
    position_robot = filtfilt(B,A,position_robot_raw);
else
    forces_mini = forces_mini_raw;
    forces_nano = forces_nano_raw;
    position_robot=position_robot_raw;
end

% Split force and acc signals
fx_mini = forces_mini(:,1);
fy_mini = forces_mini(:,2);
fz_mini = forces_mini(:,3);
tx_mini = forces_mini(:,4);
ty_mini = forces_mini(:,5);
tz_mini = forces_mini(:,6);

fx_nano = forces_nano(:,1);
fy_nano = forces_nano(:,2);
fz_nano = forces_nano(:,3);
tx_nano = forces_nano(:,4);
ty_nano = forces_nano(:,5);
tz_nano = forces_nano(:,6);

ax = acc(:,1);
ay = acc(:,2);
az = acc(:,3);

% Compute NF, TF and frame numbers during acquisition
frames = time*imAcq + 1;
nf = -(fz_mini + fz_nano);
tfy = (fy_mini + fy_nano);
tfx = (fx_mini + fx_nano);
tf = hypot(tfx,tfy);

% Pack signals in output table
signals = table(time,frames,nf,tf,tfx,tfy,...
    fx_mini,fy_mini,fz_mini,tx_mini,ty_mini,tz_mini,...
    fx_nano,fy_nano,fz_nano,tx_nano,ty_nano,tz_nano,...
    ax,ay,az,position_robot,angle_robot,cmd_robot,moisture);
end
