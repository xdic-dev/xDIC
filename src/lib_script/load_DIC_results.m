function data = load_DIC_results(varargin)
%% Parse input arguments
p = inputParser;
p.addParameter('baseDataPath',[]);
p.addParameter('baseResultPath',[]);
p.addParameter('subject',[]);
p.addParameter('material',[]);
p.addParameter('phase',[]);
p.addParameter('trial',[]);
p.addParameter('deftype',[]);

p.parse(varargin{:});
baseDataPath = p.Results.baseDataPath;
baseResultPath = p.Results.baseResultPath;
subject = p.Results.subject;
material = p.Results.material;
phase = p.Results.phase;
trial = p.Results.trial;
deftype = p.Results.deftype;

%%
fprintf('Data loading launch\n');
subject_idx = 0; 

varNames = ["subject","material","trial","phase","dir","nf","spd","spddxl","Nframe","dic"];
varTypes = ["string","string","string","string","string","double","double","string","double","struct"];
data = table('Size',[length(trial{1})*length(subject) length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

var_idx = 0; 
for subjectii = subject
    subject_idx = subject_idx+1; 
    
    % load protocol
    protocolPath = fullfile(fullfile(baseDataPath,"rawdata",subjectii,"speckles",material{subject_idx}{1},"protocol",sprintf('*.mat'))) ;
    S = dir(protocolPath);
    if isempty(S)
        error('Error: Protocol not found.');
    end
    p = load(fullfile(S.folder,S.name)); protocol = p.cond;

    dircond = protocol.table(:,strcmp(protocol.titles,'dir'));
    nfcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'nf')));
    spdcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'spd')));
    repcond = cell2mat(protocol.table(:,strcmp(protocol.titles,'rep')));
    spddxlcond = repelem(["lowspd","highspd"],length(repcond)/2); 
    ii = 0; 
    for trialii = trial{subject_idx}
        ii = ii + 1; 
        var_idx = var_idx+1; 
        trialiis = sprintf("%03.0f",trialii);
        data.subject(var_idx) = subjectii; 
        data.material(var_idx) = material{subject_idx}{ii}; 
        data.trial(var_idx) = trialiis; 
        data.phase(var_idx) = phase;
        data.dir(var_idx) = dircond{trialii};
        data.nf(var_idx) = nfcond(trialii);
        data.spd(var_idx) = spdcond(trialii);
        data.spddxl(var_idx) = spddxlcond{trialii};
        
        filename = sprintf('DIC3DPPresults_2Pairs_%s_v2.mat',deftype);
        if ~exist(fullfile(baseResultPath,'analysis',subjectii,material{subject_idx}{ii},trialiis,phase,filename),'file')
            filename = sprintf('DIC3DPPresults_2Pairs_%s_v1.mat',deftype);
        end
        
        file = load(fullfile(baseResultPath,'analysis',subjectii,material{subject_idx}{ii},trialiis,phase,filename)); 
        file = file.DIC3DPPresults; 
        
        data.Nframe(var_idx) = length(file.FaceCentroids);
        data.dic(var_idx).Faces = file.Faces; 
        data.dic(var_idx).Points3D = file.Points3D; 
        data.dic(var_idx).corrComb = file.corrComb; 
        data.dic(var_idx).FaceColors = file.FaceColors; 
        data.dic(var_idx).FacePairInds = file.FacePairInds; 
        data.dic(var_idx).PointPairInds = file.PointPairInds; 
        data.dic(var_idx).FaceCentroids = file.FaceCentroids; 
        data.dic(var_idx).FaceCorrComb = file.FaceCorrComb; 
        data.dic(var_idx).J = file.Deform.J; 
        data.dic(var_idx).Epc1 = file.Deform.Epc1; 
        data.dic(var_idx).Epc2 = file.Deform.Epc2; 
        data.dic(var_idx).EShearMax = file.Deform.EShearMax;
    end
    fprintf('done with %s\n',subjectii);
end
end