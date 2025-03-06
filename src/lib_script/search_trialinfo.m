function T = search_trialinfo(baseDataPath,subject,material_set,trial_target_set)

% Load protocol 
file_protocol = fullfile(fullfile(baseDataPath,'rawdata',subject,'speckles',material_set,...
    'protocol',sprintf('*.mat')));

S = dir(file_protocol);
p = load(fullfile(S.folder,S.name)); protocol = p.cond;

dircond = protocol.table(:,strcmp(protocol.titles,'dir'));
nfcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'nf')));
spddxlcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'spddxl')));%repelem(["lowspd","highspd"],length(repcond)/2); 

Ntrial = size(dircond,1);
Ncond = size(dircond,2);

TrialNbr = NaN(length(trial_target_set),1);
Dircond = strings(length(trial_target_set),1);
Nfcond = NaN(length(trial_target_set),1);
Spddxlcond = NaN(length(trial_target_set),1);

%correction protocol mistake
if str2double(cell2mat(regexp(subject, '\d+', 'match')))<8
    spddxlcond(1:Ntrial/2) = 0.04; 
    spddxlcond(Ntrial/2+1:Ntrial) = 0.08; 
end

for mm = 1:length(trial_target_set)
    TrialNbr(mm) = trial_target_set(mm);
    Dircond(mm) = dircond(trial_target_set(mm));
    Nfcond(mm) = nfcond(trial_target_set(mm));
    Spddxlcond(mm) = spddxlcond(trial_target_set(mm));
end
T = table(TrialNbr,Dircond,Nfcond,Spddxlcond);
end