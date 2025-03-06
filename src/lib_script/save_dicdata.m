function save_dicdata(obj, fullpathname)
% SAVE_DICDATA Saves the DIC data to a file.
%   save_dicdata(obj, fullpathname) saves the reference image, current images, and DIC data
%   from the DIC object 'obj' to the file specified by 'fullpathname'.
%
%   Inputs:
%       - obj: The DIC object containing the data to be saved.
%       - fullpathname: The full file path where the data will be saved.
%
%   Outputs:
%       None.



% Images ---------------------------------------------%
% Reference image:
reference_save = struct('type',{},'gs',{},'name',{},'path',{},'roi',{});
reference_save(1).type = obj.reference.imginfo.type;
reference_save.gs = [];
if (strcmp(obj.reference.imginfo.type,'load'))
    % Must save gs data directly
    reference_save.gs = obj.reference.imginfo.get_gs();
end
reference_save.name = obj.reference.imginfo.name;
reference_save.path = obj.reference.imginfo.path;
reference_save.roi = obj.reference.roi;

% Current image(s):
current_save = struct('type',{},'gs',{},'name',{},'path',{},'roi',{});

for i = 0:length(obj.current)-1
    current_save(i+1).type = obj.current(i+1).imginfo.type;
    current_save(i+1).gs = [];
    if (strcmp(obj.current(i+1).imginfo.type,'load'))
        % Must save gs data directly
        current_save(i+1).gs = obj.current(i+1).imginfo.get_gs();
    end
    current_save(i+1).name = obj.current(i+1).imginfo.name;
    current_save(i+1).path = obj.current(i+1).imginfo.path;
    current_save(i+1).roi = obj.current(i+1).roi;
end

% Data -----------------------------------------------%
data_dic_save = obj.data_dic;

% Save -----------------------------------------------%
save(fullpathname,'reference_save','current_save','data_dic_save','-v7.3');
end