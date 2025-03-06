function [] = plot3DmeasuresOn2Dimages(varargin)
%% fucntion for plotting 3D measures on the original 2D images
%plot3DmeasuresOn2Dimages
%plot3DmeasuresOn2Dimages(DIC2DpairResults,DIC3DpairResults)
%%
switch nargin
    case 0
        PathInitial=pwd;
%         select 3D pair result
        [fileName, filePath] = uigetfile(PathInitial,'Select the DIC3DPPresults.mat file containing the results for plotting');
        DIC3DpairResults=load([filePath fileName]);
        DIC3DPPresults=DIC3DpairResults.DIC3DPPresults;
        
        %%%%%%%%%%% select pair index %%%%%%%%%%
        answer = inputdlg('Enter pair index number');
        pairIndex=str2double(answer{1});
        
    case 2
        DIC3DPPresults=varargin{1};
        pairIndex=varargin{2};
        
    otherwise
        error('wrong number of input arguments');
end

%% select what to plot
Prompt={'\bf{Select which parameter to plot}';... % 1
    'DispMgn: Displacement magnitude (points)';... % 2
    'DispX: X Displacement (points)';... % 3
    'DispY: Y Displacement (points)';... % 4
    'DispZ:Z Displacement (points)';... % 5
    'J: Dilitation (surface area change) (faces)';... % 6
    'Lambda1: 1st principal stretch (faces)';... % 7
    'Lamda2: 2nd principal stretch (faces)';... % 8
    'Epc1: 1st (smallest) principal Lagrangian strain (faces)';... % 9
    'Epc2: 2nd (largest) principal Lagrangian strain (faces)';... % 10
    'epc1: 1st (smallest) principal Almansi strain (faces)';... % 11
    'epc2: 2nd (largest) principal Almansi strain (faces)';... % 12
    'Emgn: Lagrangian strain tensor magnitude (faces)';... % 13
    'emgn: Almansi strain tensor magnitude (faces)';... % 14
    'Eeq: Equivalent Lagrangian strain (faces)';... % 15
    'eeq: Equivalent Almansi strain (faces)';... % 16
    'EShearMax: max Lagrangian shear strain (faces)';... % 17
    'eShearMax: max Almansi shear strain (faces)';... % 18
    'Plot all';}; %19

Title='Select which parameters to plot';

Formats=struct;
Formats(1,1).type='text';
Formats(2,1).type='check';
Formats(3,1).type='check';
Formats(4,1).type='check';
Formats(5,1).type='check';
Formats(6,1).type='check';
Formats(7,1).type='check';
Formats(8,1).type='check';
Formats(9,1).type='check';
Formats(10,1).type='check';
Formats(11,1).type='check';
Formats(12,1).type='check';
Formats(13,1).type='check';
Formats(14,1).type='check';
Formats(15,1).type='check';
Formats(16,1).type='check';
Formats(17,1).type='check';
Formats(18,1).type='check';
Formats(19,1).type='check';

DefAns=cell(numel(Prompt),1);
DefAns{1}=[];
for ii=2:numel(Prompt)
    DefAns{ii}=false;
end

Options.Resize='on';
Options.FontSize=10;

[Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options);
if Answer{19}
    for ii=2:18
        Answer{ii}=true;
    end
end
if Canceled
    return
end

%% Select plotting options
optStruct=struct;
if ~Canceled && sum(cell2mat(Answer(2:end)))>0
    answer = inputdlg({'Enter maximum correlation coefficient to plot (leave blank for plotting all points/faces)'},'Input',[1,50]);
    optStruct.CorCoeffCutOff=str2double(answer{1}); % maximal correlation coefficient for display (use [] for default which is max)
else
    return
end

% %% load images
% ImPaths=DIC2DpairResults.ImPaths;
% ImSet=cell(length(ImPaths),1);
% for ii=1:length(ImPaths)
%     ImSet{ii}=imread(ImPaths{ii});
% end

%% plot animated images with the 3D measures on them
if Answer{2}
    anim8_DIC3DPP_pointMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'DispMgn',optStruct);
end
if Answer{3}
    anim8_DIC3DPP_pointMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'DispX',optStruct);
end
if Answer{4}
    anim8_DIC3DPP_pointMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'DispY',optStruct);
end
if Answer{5}
    anim8_DIC3DPP_pointMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'DispZ',optStruct);
end
if Answer{6}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'J',optStruct);
end
if Answer{7}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'Lamda1',optStruct);
end
if Answer{8}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'Lamda2',optStruct);
end
if Answer{9}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'Epc1',optStruct);
end
if Answer{10}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'Epc2',optStruct);
end
if Answer{11}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'epc1',optStruct);
end
if Answer{12}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'epc2',optStruct);
end
if Answer{13}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'Emgn',optStruct);
end
if Answer{14}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'emgn',optStruct);
end
if Answer{15}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'Eeq',optStruct);
end
if Answer{16}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'eeq',optStruct);
end
if Answer{17}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'EShearMax',optStruct);
end
if Answer{18}
    anim8_DIC3DPP_faceMeasure_onImages_n_n(DIC3DPPresults,pairIndex,'eShearMax',optStruct);
end


end

 
%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>