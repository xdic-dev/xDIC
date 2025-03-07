%% Add MultiDIC library path so functions are known to use here

mPath = pwd;
addpath(genpath(fullfile(mPath,'toolbox','MultiDIC-master','lib_MultiDIC'))); % add libraries
addpath(genpath(fullfile(mPath,'toolbox','MultiDIC-master','lib_ext'))); % add external library
addpath(genpath(fullfile(mPath,'toolbox','MultiDIC-master','main_scripts'))); % add main script 
addpath(genpath(fullfile(mPath,'toolbox','updateABRMultiDIC'))); % add rewrited code from Arbaud Browet
savepath
% install Ncorr
cd(fullfile(mPath,'toolbox','MultiDIC-master','lib_ext','ncorr_2D_matlab-master'));
handles_ncorr=ncorr;
cd(mPath);

% message
h=msgbox('MultiDIC installed successfully');
hp=get(h, 'position');
set(h, 'position', [hp(1) hp(2) 150 50]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set(ch, 'FontSize',10);