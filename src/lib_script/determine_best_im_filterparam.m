
baseDataPath = "E:\transfert_SMP"; 
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 
addpath(genpath(fullfile('script'))); % add external library

% Frictional conditions
FRcond{1} = "glass"; 
FRcond{2} = "coating"; 

%global parameters
Npair = 2; %Number of pair 
Fs = 1e3; %robot sampling frequency
fps = 50; %image sampling frequency 

subject = "S01"; 
phase = "loading";
material_nbr = 2; 
showvisu = 1; %boolean for visualization

trialnum = 0; 
param_init = [20 150];
param_filt_im = [param_init+25*0;
                 param_init+25*1;
                 param_init+25*2;
                 param_init+25*3;
                 param_init+25*4 ];

for jj = 1:length(param_filt_im)
    param_filt_im_jj = param_filt_im(jj,:); 
    trialnum = trialnum+1;
    % frame index for analysis
    if strcmp(phase,'loading')
        frameidx = load(fullfile('analysis',subject,'frameidx.mat')); frameidx = frameidx.frameidx;
        idx_frame_start = frameidx{material_nbr}.(phase).start(1);
        idx_frame_end = frameidx{material_nbr}.(phase).end(1);
    else
        idx_frame_start = [];
        idx_frame_end = [];
    end
    calibPath = fullfile(baseResultPath,"analysis",subject,"calib","4");
    
    [outputPath,pairOrder,pairForced] = stepD_2DDIC(...
        'baseDataPath',baseDataPath,...
        'baseResultPath',baseResultPath,...
        'subject',subject,...
        'material',FRcond{material_nbr},...
        'trial',sprintf("%03d",1),...
        'stereopair',1,...
        'phase',phase,...
        'showvisu',showvisu,...
        'idxstart_set',idx_frame_start,...
        'idxend_set',idx_frame_end,...
        'param_filt_im',param_filt_im_jj);
    stepE_3DReconstruction('basePath',outputPath,...
        'calibPath',calibPath,...
        'showvisu',showvisu,...
        'pairOrder',pairOrder,...
        'pairForced',pairForced,...
        'title_ax',sprintf("%s-%s-%s-%s",subject,FRcond{material_nbr},sprintf("%03d",trial_jj),phase));

end