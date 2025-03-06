function trials = search_trial2target(baseDataPath,subject,phase,material_set,nfcond_set,spddxlcond_set)

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
%correction protocol mistake
if str2double(cell2mat(regexp(subject, '\d+', 'match')))<8
    spddxlcond(1:Ntrial/2) = 0.04; 
    spddxlcond(Ntrial/2+1:Ntrial) = 0.08; 
end

trials = []; 
trialnum = 1:Ntrial; 
for ii = 1:length(nfcond_set)
    nfcond_set_ii = nfcond_set(ii); 
    for jj = 1:length(spddxlcond_set)
        spddxlcond_set_jj = spddxlcond_set(jj);
        if strcmp(phase,'loading') 
            trialLogic = (strcmp(dircond,'Ubnf')|strcmp(dircond,'Rbnf'))...
                &(nfcond==nfcond_set_ii&spddxlcond==spddxlcond_set_jj); 
        else 
            trialLogic = true(Ntrial,Ncond); 
        end
        trials = cat(2,trials,trialnum(trialLogic));
    end
end
end