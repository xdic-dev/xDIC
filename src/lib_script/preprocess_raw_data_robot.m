%%% PREPROCESS RAW DATA ROBOT
%% Global parameters
subject = 'S05';

%%
FRcond = {"glass","coating"};

color_FRcond_NF{1,1} = rgb('Salmon');
color_FRcond_NF{2,1} = rgb('OliveDrab');
color_FRcond_NF{1,2} = rgb('Salmon');
color_FRcond_NF{2,2} = rgb('OliveDrab');

Fs = 1e3; 
Filterfrequency = 10; 
fps = 50; 
ref_coef = 0.05;
baseline_threshold_touch = 0.5;
range_CF_eval = 0.2; % fraction of the tangential phase  
nbrimagemargin = 5;
Nframe = 75; 

%
width_figure = 20;

%
num = regexp(subject, '\d+', 'match');
subject_num = str2double(num{end});

go2ref_phase_number = 1; 
waiting_phase_number = 2; 
loading_phase_number = 3; 
sliding1_phase_number = 4; 
relax1_phase_number = 5; 
sliding2_phase_number = 6;
relax2_phase_number = 7; 

%% Import
data_robot = struct(); 
protocol = struct(); 
for pp = 1:length(FRcond) %FR cond
    
    %protocol file loading 
    protocolPath = fullfile(baseDataPath,"rawdata",subject,fingerpattern,FRcond{pp},"protocol",sprintf('*.mat')) ;
    S = dir(protocolPath);
    if isempty(S)
        error('Error: Protocol not found.');
    else 
        disp('--> protocol loaded');
    end
    p = load(fullfile(S.folder,S.name)); protocol.(FRcond{pp}) = p.cond;
    Ntrial = size(protocol.(FRcond{pp}).table,1);
    
    %raw data robot files loading
    datafilepath = fullfile(baseDataPath,"rawdata",subject,fingerpattern,FRcond{pp},"robot");
    [data_robot.(FRcond{pp}),missing_files] = ...
        import_robotdata_subject(datafilepath,'FreqFilter',Filterfrequency,'Ntrial',Ntrial); 
    
    if missing_files
        p = protocol.(FRcond{pp}); 
        p_spdcondmask = cell2mat(p.table(:,strcmp(p.titles,'spd')));
        emptyCells = cellfun(@isempty,data_robot.(FRcond{pp}));
        idx_missing_data = find(emptyCells); 
        for idx = idx_missing_data
            spdcond_idx = p_spdcondmask(idx); 
            %take the first iteration with the time course
            d = data_robot.(FRcond{pp}){find(p_spdcondmask == spdcond_idx,1)}; 
            %fill NaN table with the same structure
            T = array2table(NaN(size(d)-[0,4])); %minus table size which are emebedded inside the large table
            T.Properties.VariableNames = d.Properties.VariableNames(1:end-4);
            T = [T,table(NaN(size(d,1),3),NaN(size(d,1),3),NaN(size(d,1),3),NaN(size(d,1),1),...
                'VariableNames',d.Properties.VariableNames(end-3:end))];  
            data_robot.(FRcond{pp}){idx} = T; 
        end
    end
end


%% Protocol information : Index formatting and timings
pmask = cell(length(FRcond),1); 
for pp = 1:length(FRcond) %FR cond
    p = protocol.(FRcond{pp}); 
    pmask{pp}.dircond_mask = convertCharsToStrings(p.table(:,strcmp(p.titles,'dir')));
    pmask{pp}.spdcond_mask = cell2mat(p.table(:,strcmp(p.titles,'spd')));
    pmask{pp}.nfcond_mask = cell2mat(p.table(:,strcmp(p.titles,'nf')));
    pmask{pp}.repcond_mask = cell2mat(p.table(:,strcmp(p.titles,'rep')));
    pmask{pp}.spddxlcond_mask = cell2mat(p.table(:,strcmp(p.titles,'spddxl')));
    
    pmask{pp}.dircond = unique(pmask{pp}.dircond_mask)'; 
    pmask{pp}.spdcond = unique(pmask{pp}.spdcond_mask)'; 
    pmask{pp}.nfcond = unique(pmask{pp}.nfcond_mask)';
    pmask{pp}.repcond = unique(pmask{pp}.repcond_mask)'; 
    pmask{pp}.spddxlcond = unique(pmask{pp}.spddxlcond_mask)';
    
    if subject_num < 8
        %correction
        pmask{pp}.spddxlcond_mask = repelem(pmask{pp}.spddxlcond,length(p.table)/length(pmask{pp}.spddxlcond))'; %error in the protocol
    end
    % Time scale for the friction coefficient estimation
    timing_protocol = cumsum(p.dur);
    pmask{pp}.loading_time = cell(length(pmask{pp}.spddxlcond),1);
    pmask{pp}.sliding1_time = cell(length(pmask{pp}.spdcond),1);   
    pmask{pp}.relax1_time = cell(length(pmask{pp}.spdcond),1); 
    pmask{pp}.sliding2_time = cell(length(pmask{pp}.spdcond),1); 
    pmask{pp}.relax2_time = cell(length(pmask{pp}.spdcond),1); 
    for jj = 1:length(pmask{pp}.spdcond)
        pmask{pp}.loading_time{jj} = timing_protocol(loading_phase_number-1)+1:timing_protocol(loading_phase_number);
        pmask{pp}.sliding1_time{jj} = timing_protocol(sliding1_phase_number-1)+1+zeros(length(pmask{pp}.spdcond_mask),1):...
            timing_protocol(sliding1_phase_number-1)+p.dst/pmask{pp}.spdcond(jj)*Fs;
        pmask{pp}.relax1_time{jj} = pmask{pp}.sliding1_time{jj}(end)+1:...
            pmask{pp}.sliding1_time{jj}(end)+timing_protocol(relax1_phase_number)-timing_protocol(relax1_phase_number-1);
        pmask{pp}.sliding2_time{jj} = pmask{pp}.relax1_time{jj}(end)+1:...
            pmask{pp}.relax1_time{jj}(end)+p.dst/pmask{pp}.spdcond(jj)*Fs;
        pmask{pp}.relax2_time{jj} = pmask{pp}.sliding2_time{jj}+1:...
            pmask{pp}.sliding2_time{jj}(end)+timing_protocol(relax1_phase_number)-timing_protocol(relax1_phase_number-1); 
    end
end
% % Finding touch instant 
% idx = cell(1,length(FRcond));
% for pp = 1:length(FRcond) %FR cond
%     idx{pp}.loading_start = NaN(1,length(p.table)); 
%     idx{pp}.loading_end = NaN(1,length(p.table)); 
%     for itrial = 1:length(pmask{pp}.dircond_mask) % over all trials
%         spdcond_itrial = pmask{pp}.spdcond_mask(itrial)==pmask{pp}.spdcond; 
%         nfcond_itrial = pmask{pp}.nfcond_mask(itrial)==pmask{pp}.nfcond; 
%         
%         nf = data_robot.(FRcond{pp}){itrial}.nf; 
%         
%         idx_load_timestart = pmask{pp}.loading_time{spdcond_itrial}(1);
%         idx_load_timeend = pmask{pp}.loading_time{spdcond_itrial}(end);
%         
%         idx_estim_start = idx_load_timestart+find(nf(idx_load_timestart:idx_load_timeend)>=baseline_threshold_touch,1);
%         mean_baseline = mean(nf(idx_load_timestart:idx_estim_start-100));
%         
%         idx{pp}.loading.start(itrial) = find(nf(idx_load_timestart:idx_load_timeend)>=mean_baseline+baseline_threshold_touch,1);
%         idx{pp}.loading.end(itrial) = find(nf(idx_load_timestart:idx_load_timeend)>=nfcond_itrial*(1-ref_coef),1);
%         
%     end
% end
%% NF offset
% We take the average data over 500ms before the loading phase.
figure; 
ax1 = subplot(4,1,1); hold on; 
ax2 = subplot(4,1,2); hold on; 
ax3 = subplot(4,1,3); hold on; 
ax4 = subplot(4,1,4); hold on; 
Ntrial = length(data_robot.(FRcond{1})); 
Nspdcond = length(pmask{pp}.spdcond);
Nfrcond = length(FRcond);
nf_all = cell(Nfrcond,Nspdcond);

invalid_trials_idx = []; 
for pp = 1:Nfrcond %FR cond
for itrial = 1:Ntrial 
idxtimestart = pmask{pp}.loading_time{1}(1); 
temp1 = data_robot.(FRcond{pp}){itrial}.time;
temp2 = data_robot.(FRcond{pp}){itrial}.nf;
baseline1_mean = mean(temp2(idxtimestart-500:idxtimestart,:),1);
temp2 = temp2 - baseline1_mean; 

%check whether something was touching the plate just before loading 
if abs(baseline1_mean) > 0.5
    invalid_trials_idx = itrial+Ntrial*(pp-1); 
end
%find the start of actual loading 
idx = find(temp2(idxtimestart:end)>baseline_threshold_touch,1)+idxtimestart; 
idx = round(idx/Fs*fps)-nbrimagemargin;
idx = idx/fps*Fs; 

%plot
plot(ax1,temp1,data_robot.(FRcond{pp}){itrial}.nf);
plot(ax2,temp1,temp2); 
plot(ax2,temp1(idx), temp2(idx),'rx');
plot(ax3,temp1(idx:end)-temp1(idx), temp2(idx:end)); 
plot(ax4,temp1(idx:(idx+Nframe/fps*Fs))-temp1(idx), temp2(idx:(idx+Nframe/fps*Fs))); 


% %output var
% data_robot.(FRcond{pp}){itrial}.nf = temp2; 

end
end
fillgray(ax2,temp1(idxtimestart-500),temp1(idxtimestart))
xline(ax2,[temp1(idxtimestart-500),temp1(idxtimestart)]);
yline(ax2,baseline_threshold_touch);
hold off; 


%% Wrong trials 

for itrial = invalid_trials_idx
    if itrial>Ntrial
        itrial = itrial-Ntrial; 
        pp = 2; 
    else 
        pp = 1; 
    end
    p = protocol.(FRcond{pp});
    p_spdcondmask = cell2mat(p.table(:,strcmp(p.titles,'spd')));
    
    spdcond_idx = p_spdcondmask(itrial);
    %take the first iteration with the time course
    d = data_robot.(FRcond{pp}){find(p_spdcondmask == spdcond_idx,1)};
    %fill NaN table with the same structure
    T = array2table(NaN(size(d)-[0,4])); %minus table size which are emebedded inside the large table
    T.Properties.VariableNames = d.Properties.VariableNames(1:end-4);
    T = [T,table(NaN(size(d,1),3),NaN(size(d,1),3),NaN(size(d,1),3),NaN(size(d,1),1),...
        'VariableNames',d.Properties.VariableNames(end-3:end))];
    data_robot.(FRcond{pp}){itrial} = T;
end

%%
% %make directory 
% folder_loc = fullfile(baseResultPath,"robot_process");
% if ~exist(folder_loc,'dir')
%     mkdir(folder_loc);
% end
% filename = 'robot_'+subject;
% file_loc = fullfile(folder_loc,filename+'.mat'); 
% save(file_loc+'.mat')

