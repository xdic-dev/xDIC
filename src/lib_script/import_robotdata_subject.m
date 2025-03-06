function varargout = import_robotdata_subject(basePath,varargin)
%% input 
% Parse input arguments
p = inputParser;
p.addRequired('basePath');
p.addParameter('FreqFilter',40);
p.addParameter('filenumber',[]);
p.addParameter('Ntrial',[]);

% Parse input arguments
p.parse(basePath,varargin{:});
FreqFilter = p.Results.FreqFilter;
filenumber= p.Results.filenumber;
Ntrial= p.Results.Ntrial;

%% 
%directory = "rawdata/S01/robot/speckles";
directory = fullfile(basePath);
filesAndFolders = dir(directory);     % Returns all the files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));

numOfFiles = length(filesInDir);
if isempty(filenumber)
    filenumber=1:numOfFiles;
end

if isempty(Ntrial)
    Ntrial = numOfFiles; 
end

% Warning statements
missing_files = false; 
if (Ntrial > numOfFiles)
    disp('Warning : The number of files found does not match the protocol number of trials'); 
    disp('          NaN will be filled for missing trials');
    missing_files = true; 
elseif (Ntrial < numOfFiles)
    error('Warning : The number of files found does not match the protocol number of trials\n Clean your folder!');
end

fprintf("\n---> %d files to be processed : \n",numOfFiles);

data = cell(Ntrial,1);
for ii = filenumber
    
    filename = filesInDir(ii).name;
    S = regexp(filename,'_','split');
    if ii == 1
        fprintf("protocol : "); 
        fprintf("%s ",S{1:end-4})
        fprintf("\n"); 
        fprintf("%s",repelem("|",numOfFiles)); 
        fprintf("\n"); 
    end
    num = regexp(filename, '\d+', 'match');
    idx = str2double(num{end});
%     file = sprintf("%s ",S{end-3:end}); 
%     fprintf("file %s in proccess num %d\n",filename,idx);
    fprintf("|"); 
    
    %import table data into cells corresponding to the right trial
    data{idx} = smp_import(fullfile(directory,filename),1,FreqFilter);
end
    
fprintf("\n");
fprintf("---> data loaded\n");

%% output 
varargout{1} = data;

if nargout>1 
    varargout{2} = missing_files; 
end

end