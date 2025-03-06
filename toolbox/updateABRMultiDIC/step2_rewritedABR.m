%% parameters

% output path
outputPath='/Users/abrowet/finger_dic/acquisition_march/testNcorrMerge/results/sub_trial008';

% path to video frames
% basePath = "/Users/abrowet/finger_dic/acquisition_march/raw/speckles/frames/DD_dic_006_525_050_trial008_";
basePath = "/Users/abrowet/finger_dic/acquisition_march/raw/speckles/frames/sub_trial008";

% ensure output path exists
if (~ exist(outputPath, "dir"))
    disp('creating output path')
    mkdir(outputPath)
end

% index and path to 1st camera
cam_1 = 1;
cam_first = sprintf("cam_%d", cam_1);
cam_first_path = fullfile(basePath, cam_first);

% index and path to 2nd camera
cam_2 = 4;
cam_second = sprintf("cam_%d", cam_2);
cam_second_path = fullfile(basePath, cam_second);

% image extension
image_ext = '.jpg';

% Ncorr DIC parameters
step1_parameters = struct( ...
    'type', 'regular',...
    'radius', 45,...
    'spacing', 10, ...
    'cutoff_diffnorm', 1e-6, ...
    'cutoff_iteration', 100, ...
    'total_threads', 1, ...
    'stepanalysis_params', struct( ...
                            'enabled', 1, ...
                            'type', 'seed', ...
                            'auto',1, ...
                            'step', 5) ...
    );

step1_2_parameters = struct( ...
    'type', 'regular',...
    'radius', 45,...
    'spacing', 10, ...
    'cutoff_diffnorm', 1e-6, ...
    'cutoff_iteration', 100, ...
    'total_threads', 1, ...
    'stepanalysis_params', struct( ...
                            'enabled', 0, ...
                            'type', 'seed', ...
                            'auto',1, ...
                            'step', 5) ...
    );

step2_parameters = step1_parameters;

% specific parameters
% !!! initial_seed parameter is set according to the original image size
step1_parameters.initial_seed = [1050, 520];

disp('Done')
%% Check folder content
% find the number of images
cam1_files = dir(fullfile(cam_first_path,sprintf('/*%s',image_ext)));
nb_images = size(cam1_files,1);

cam2_files = dir(fullfile(cam_second_path,sprintf('/*%s',image_ext)));
nb_images_2 = size(cam2_files,1);

if (nb_images ~= nb_images_2)
    error('Your folders contain different number of images')
end

disp('Done')
% nb_images = 2;
% nb_images_2 = 2;

%% read 1st set
disp('reading camera 1 images')
set1 = cell(nb_images,1);
for i = 1:nb_images
    % disp(i)
    ii = imread( ...
        fullfile(basePath, cam_first , ...
        strjoin(["frame_" sprintf('%04d',i) ".jpg"],"")) ...
        );
    if size(ii,3)==3
        ii = Irgb2gray(ii);
    end
    set1{i} =  uint8(ii);
end
%set1=uint8(set1);
disp('  done')

%%
if (length(step1_parameters.initial_seed)==2)
    close all
    imshow(set1{1})
    hold on
    plot(step1_parameters.initial_seed(1), ...
        step1_parameters.initial_seed(2), ...
        '*r')
end

%% start ncorr to create ROI
handles = ncorr;
handles.set_ref(set1{1});
handles.set_cur(set1{1});
disp('Draw ROI ')

%% set ROI file path
ROIfile = fullfile( outputPath, "ROI.mat");

%% save ROI 
ROImask = handles.reference.roi;
save(ROIfile, "ROImask", '-v7.3');

close(handles.handles_gui.figure);
disp('save ROI done');


%% run ncorr on 1st
showGui = false;
if (isfield(step1_parameters, "initial_seed"))
    initial_seed = step1_parameters.initial_seed;
else
     initial_seed = [];
end
handles = ncorr_abr(showGui, initial_seed);

% handles = ncorr;
handles.set_ref(set1{1});
handles.set_cur(set1);
ROI = load(ROIfile);
handles.set_roi_ref(ROI.ROImask.mask); 

%%
if (true)
    handles.set_dic_parameters_manually( ...
        step1_parameters.type, ...
        step1_parameters.radius, ...
        step1_parameters.spacing, ...
        step1_parameters.cutoff_diffnorm, ...
        step1_parameters.cutoff_iteration,...
        step1_parameters.total_threads, ...
        step1_parameters.stepanalysis_params)

    % refresh the handle
    handles.abr_refresh()
    
    % trigger DIC
    handles.callback_topmenu_dic();

    % trigger Displacement
    handles.callback_topmenu_formatdisp();

end

disp('Done')
% SET ROI IN NCORR => DONE
% SET DIC PARAMS (WITH Seed Propagation) => DONE
% PERFORM DIC => DONE
% FORMAT DISPLACEMENT => DONE

%% save ncorr results
file1 = fullfile( outputPath, "ncorr1.mat");
data_dic1 = handles.data_dic;

% that's a bit too much to save
data_dic1.reference = handles.reference;
% data_dic1.current = handles.current;


save(file1, "data_dic1", '-v7.3');

close(handles.handles_gui.figure);

ref1_roi = data_dic1.reference.roi;
clear handles
disp('save done');



%% load only the 1st image on the second set
set2 = set1(1);
ii = imread( ...
    fullfile(basePath, cam_second ,"frame_0001.jpg") ...
    );
if size(ii,3)==3
    ii = Irgb2gray(ii);
end
set2{2} =  uint8(ii);
disp('reading 1st image of set 2 done');

%% run ncorr on the firsts images of each set 
% WE KEEP THE ROI OF THE FIRST SET
showGui = true;
handles = ncorr_abr(showGui);
% handles = ncorr();
handles.set_ref(set2{1});
handles.set_cur(set2(2));

handles.set_roi_ref(ref1_roi.mask);

if (true)
    handles.set_dic_parameters_manually( ...
        step1_2_parameters.type, ...
        step1_2_parameters.radius, ...
        step1_2_parameters.spacing, ...
        step1_2_parameters.cutoff_diffnorm, ...
        step1_2_parameters.cutoff_iteration,...
        step1_2_parameters.total_threads, ...
        step1_2_parameters.stepanalysis_params)
    % refresh the handle
    handles.abr_refresh()
    
    % trigger DIC
    handles.callback_topmenu_dic();

    % trigger Displacement
    handles.callback_topmenu_formatdisp();

end
disp('Done')

% SET DIC PARAMS (NO Seed Propagation)
% PERFORM DIC
% FORMAT DISPLACEMENT

%% save ncorr results
file1_2 = fullfile( outputPath, "ncorr1_2.mat");
data_dic1_2 = handles.data_dic;

data_dic1_2.current = handles.current;

ref2_roi = handles.current.roi;
save(file1_2, "data_dic1_2");

close(handles.handles_gui.figure);
disp('save done');

%% show the2 roi
figure
subplot(2,2,1)
imshow(set2{1})
subplot(2,2,2)
imshow(set2{2})

subplot(2,2,3)
mm = set2{1}.*uint8(ref1_roi.mask);
mm(mm==0) = 100;
imshow(mm)
subplot(2,2,4)
mm=set2{2}.*uint8(ref2_roi.mask);
mm(mm==0) = 100;
imshow(mm);



%% read 2nd set
set2 = cell(nb_images,1);
for i = 1:nb_images
    % disp(i)
    ii = imread( ...
        fullfile(basePath, cam_second , ...
        strjoin(["frame_" sprintf('%04d',i) ".jpg"],"")) ...
        );
    if size(ii,3)==3
        ii = Irgb2gray(ii);
    end
    set2{i} =  uint8(ii);
end
disp('reading done')



%% run ncorr on the second set of images
% WE KEEP THE ROI OF THE FIRST SET
% handles = ncorr;
% handles.set_ref(set2{1});
% handles.set_cur( set2(2:end));
% %roi_mask =
% handles.set_roi_ref(ref2_roi.mask);


% WE KEEP THE ROI OF THE FIRST SET
showGui = false;
handles = ncorr_abr(showGui);
% handles = ncorr;
handles.set_ref(set2{1});
handles.set_cur( set2(2:end));

handles.set_roi_ref(ref2_roi.mask);

if (true)
    handles.set_dic_parameters_manually( ...
        step2_parameters.type, ...
        step2_parameters.radius, ...
        step2_parameters.spacing, ...
        step2_parameters.cutoff_diffnorm, ...
        step2_parameters.cutoff_iteration,...
        step2_parameters.total_threads, ...
        step2_parameters.stepanalysis_params)
    
    % refresh the handle
    handles.abr_refresh()

    % trigger DIC
    handles.callback_topmenu_dic();

    % trigger Displacement
    handles.callback_topmenu_formatdisp();

end

disp('Done')

% SET DIC PARAMS (WITH Seed Propagation)
% PERFORM DIC
% FORMAT DISPLACEMENT


%% save ncorr results
file2 = fullfile( outputPath, "ncorr2.mat");
data_dic2 = handles.data_dic;

% data_dic1_2.current = handles.current;
%ref2_roi = handles.current.roi;

save(file2, "data_dic2", '-v7.3');

close(handles.handles_gui.figure);
disp('save done');


%% clear all then perform DIC
%clear;

ncorr1 = load(fullfile( outputPath, "ncorr1.mat"));
data1 = ncorr1.data_dic1;

ncorr1_2 = load(fullfile( outputPath, "ncorr1_2.mat"));
data1_2 = ncorr1_2.data_dic1_2;

ncorr2 = load(fullfile( outputPath, "ncorr2.mat"));
data2 = ncorr2.data_dic2;

disp('loading Done')

%% start DIC data structure
DIC_results= struct;

DIC_results.nCamRef = cam_1;
DIC_results.nCamDef = cam_2;
DIC_results.nImages = size(data1.displacements,2);

DIC_results.ImPaths = {};
for cam = [DIC_results.nCamRef, DIC_results.nCamDef]
    for im = 1:DIC_results.nImages
        o = sprintf('%scam_%d/frame_%04d.jpg',basePath, cam, im);
        DIC_results.ImPaths{end+1,1} = o;
    end
end

DIC_results.ROImask = data1.reference.roi.mask;

DIC_results.ncorrInfo = data1.dispinfo;
DIC_results.ncorrInfo.cutoff_corrcoef = [
    DIC_results.ncorrInfo.cutoff_corrcoef;
    data1_2.dispinfo.cutoff_corrcoef;
    data2.dispinfo.cutoff_corrcoef;
    ];

disp('init Done')

%% compute points / corCoeffVec / Faces / FaceColors

% Extract results
Disp1=data1.displacements;
DispInfo1=data1.dispinfo;

Disp1_2=data1_2.displacements;
DispInfo1_2=data1_2.dispinfo;

Disp2=data2.displacements;
DispInfo2=data2.dispinfo;

nCur=size(Disp1,2) + 1 + size(Disp2,2);

Factor=DispInfo1.spacing+1;

% Define output structure
ROI_DIC=cell(nCur,1);
CorCoeff=cell(nCur,1); 
CorCoeffVec=cell(nCur,1); 
Points=cell(nCur,1);
Uvec=cell(nCur,1); 
Vvec=cell(nCur,1);


[YrefROIVec1,XrefROIVec1] = find(Disp1(1).roi_dic.mask);
PtempRef1=[XrefROIVec1,YrefROIVec1];
PtempRef1=(PtempRef1-1)*Factor+1; % switch from sapcing to pixels
Pref=PtempRef1;

IMref = imread( ...
        fullfile(basePath, cam_first , ...
        strjoin(["frame_" sprintf('%04d',1) ".jpg"],"")) ...
        );
if size(IMref,3)==3
    IMref = Irgb2gray(IMref);
end

disp('data ready')

%% First loop on data1
outside_it = 0;
for ii=1:size(Disp1,2)

    outside_it = outside_it+1;
    
    ROI_DIC{ii}=Disp1(ii).roi_dic.mask;
    
    CorCoeff{ii}=Disp1(ii).plot_corrcoef_dic;
    Uref=Disp1(ii).plot_u_ref_formatted;
    Vref=Disp1(ii).plot_v_ref_formatted;
    
%     [YrefROIVec,XrefROIVec] = find(ROI_DIC{1});
%     [YcurROIVec,XcurROIVec] = find(ROI_DIC{ii});
%     
%     PtempRef=[XrefROIVec,YrefROIVec];
%     PtempRef=(PtempRef-1)*Factor+1; % switch from sapcing to pixels
%     Pref=PtempRef;
       
    CorCoeffVec{ii}=CorCoeff{ii}(ROI_DIC{1});
    CorCoeffVec{ii}(CorCoeffVec{ii}==0)=NaN;
    
    % displacements from ref to cur
    UrefROIVec=Uref(ROI_DIC{1});
    UrefROIVec(UrefROIVec==0)=NaN;
    VrefROIVec=Vref(ROI_DIC{1});
    VrefROIVec(VrefROIVec==0)=NaN;
    
    Uvec{ii}=UrefROIVec;
    Vvec{ii}=UrefROIVec;
    % current points
    Points{ii}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];
    
    % save face colors for further 3D analysis
    if ii==1
        % pixel colors
        IMrefSmall=IMref(1:Factor:end,1:Factor:end);
        IMrefSmallMasked=IMrefSmall;
        IMrefSmallMasked(~ROI_DIC{1})=[];
        ColorRef=IMrefSmallMasked(:);
    end
    
end

disp('First loop done')


%% Match between ref of the 2 sets

outside_it=outside_it+1;
ii=outside_it;

ROI_DIC{ii}=Disp1_2(1).roi_dic.mask;

CorCoeff{ii}=Disp1_2(1).plot_corrcoef_dic;
Uref=Disp1_2(1).plot_u_ref_formatted;
Vref=Disp1_2(1).plot_v_ref_formatted;
   
U1_2 = Uref;
V1_2 = Vref;

CorCoeffVec{ii}=CorCoeff{ii}(ROI_DIC{1});
CorCoeffVec{ii}(CorCoeffVec{ii}==0)=NaN;

% displacements from ref to cur
UrefROIVec=Uref(ROI_DIC{1});
UrefROIVec(UrefROIVec==0)=NaN;
VrefROIVec=Vref(ROI_DIC{1});
VrefROIVec(VrefROIVec==0)=NaN;

Uvec{ii}=UrefROIVec;
Vvec{ii}=UrefROIVec;
% current points
Points{ii}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];
disp('Matching 1-2 done')

%% compute point mapping
U_mapped = U1_2/Factor;
V_mapped = V1_2/Factor;

[U_mapped_row, U_mapped_col, U_mapped_value] = find(U_mapped);
[V_mapped_row, V_mapped_col, V_mapped_value] = find(V_mapped);

mapped_indices = sub2ind(size(U_mapped), U_mapped_row, U_mapped_col);

U_down = floor(U_mapped_value);
U_up = ceil(U_mapped_value);
V_down = floor(V_mapped_value);
V_up = ceil(V_mapped_value);

U_low = U_mapped_col + U_down;
U_high = U_mapped_col + U_up;
V_low = V_mapped_row + V_down;
V_high = V_mapped_row + V_up;

ind_low_low = sub2ind(size(U_mapped), V_low, U_low);
ind_high_low = sub2ind(size(U_mapped), V_low, U_high);
ind_low_high = sub2ind(size(U_mapped), V_high, U_low);
ind_high_high = sub2ind(size(U_mapped), V_high, U_high);

real_col = U_mapped_col + U_mapped_value;
real_row = V_mapped_row + V_mapped_value;

% compute distances
s(:,1) = 1./sqrt((real_col-U_low).^2 + (real_row-V_low).^2);
s(:,2) = 1./sqrt((real_col-U_high).^2 + (real_row-V_low).^2);
s(:,3) = 1./sqrt((real_col-U_low).^2 + (real_row-V_high).^2);
s(:,4) = 1./sqrt((real_col-U_high).^2 + (real_row-V_high).^2);


%% Second loop on data2
for ii=1:size(Disp2,2)
    
    outside_it = outside_it+1;
    ROI_DIC{outside_it}=Disp2(ii).roi_dic.mask;
    
    CorCoeff{outside_it}=Disp2(ii).plot_corrcoef_dic;
    Uref=Disp2(ii).plot_u_ref_formatted;
    Vref=Disp2(ii).plot_v_ref_formatted;

    all_u = [
        Uref(ind_low_low) Uref(ind_high_low) Uref(ind_low_high) Uref(ind_high_high)
        ];
    %all_u(all_u==0)=nan;

    all_v = [Vref(ind_low_low) Vref(ind_high_low) Vref(ind_low_high) Vref(ind_high_high)
        ];
    %all_v(all_v==0)=nan;

    scaled_u = sum(all_u.*s, 2)./sum((all_u~=0).*s, 2);
    scaled_v = sum(all_v.*s, 2)./sum((all_v~=0).*s, 2);
    
    scaled_u(scaled_u==0)=nan;
    scaled_v(scaled_v==0)=nan;

    new_Uref = zeros(size(Uref));
    new_Vref = zeros(size(Vref));
    
    new_Uref(mapped_indices) = scaled_u;
    new_Vref(mapped_indices) = scaled_v;

%     p_row = 35;
%     p_col=71;
% 
%     shift_x = U_mapped(p_row, p_col);
%     shift_y = V_mapped(p_row, p_col);
% 
%     new_p_row = p_row + shift_y;
%     new_p_col = p_col + shift_x;
% 
%     u_mean = mean([
%     Uref(floor(new_p_row), floor(new_p_col))
%     Uref(ceil(new_p_row), floor(new_p_col))
%     Uref(floor(new_p_row), ceil(new_p_col))
%     Uref(ceil(new_p_row), ceil(new_p_col))
%     ]);
%     v_mean = mean([
%     Vref(floor(new_p_row), floor(new_p_col))
%     Vref(ceil(new_p_row), floor(new_p_col))
%     Vref(floor(new_p_row), ceil(new_p_col))
%     Vref(ceil(new_p_row), ceil(new_p_col))
%     ]);
%     disp([Uref(p_row, p_col) u_mean])
%     disp([Vref(p_row, p_col) v_mean])
%     disp( [new_Uref(p_row, p_col) new_Vref(p_row, p_col)])

    % correct the displacement
    Uref = new_Uref + U1_2;
    Vref = new_Vref + V1_2;
    
%     [YrefROIVec,XrefROIVec] = find(ROI_DIC{1});
%     [YcurROIVec,XcurROIVec] = find(ROI_DIC{ii});
%     
%     PtempRef=[XrefROIVec,YrefROIVec];
%     PtempRef=(PtempRef-1)*Factor+1; % switch from sapcing to pixels
%     Pref=PtempRef;
       
    CorCoeffVec{outside_it}=CorCoeff{outside_it}(ROI_DIC{1});
    CorCoeffVec{outside_it}(CorCoeffVec{outside_it}==0)=NaN;
    
    % displacements from ref to cur
    UrefROIVec=Uref(ROI_DIC{1});
    UrefROIVec(UrefROIVec==0)=NaN;
    VrefROIVec=Vref(ROI_DIC{1});
    VrefROIVec(VrefROIVec==0)=NaN;
    
    Uvec{outside_it}=UrefROIVec;
    Vvec{outside_it}=UrefROIVec;
    % current points
    Points{outside_it}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];    
end

disp('Second loop done')

%% Finish step 2
% Create triangulation
DT = delaunayTriangulation(Pref);
F = DT.ConnectivityList;
V = Pref;

% remove irregular triangles
EdgeLengths = patchEdgeLengths(F,V);
EdgeLengths = [EdgeLengths(1:3:length(EdgeLengths)) EdgeLengths(2:3:length(EdgeLengths)) EdgeLengths(3:3:length(EdgeLengths))];
EdgeLengthsMax=max(EdgeLengths,[],2);

F(EdgeLengthsMax>1.1*sqrt(2)*Factor,:)=[];

% flip direction to have normals pointing out
F=F(:,[1 3 2]);

% face colors (average node colors)
CF=mean(ColorRef(F),2);

disp('Triangulation Done');

%% Add data to output structure & save
DIC_results.Points=Points;
DIC_results.CorCoeffVec=CorCoeffVec;
DIC_results.Faces=F;
DIC_results.FaceColors=CF;


output_filename = fullfile( outputPath, ...
    sprintf("myDIC2DpairResults_C_%d_C_%d.mat", cam_1, cam_2));

DIC2DpairResults = DIC_results;
save(output_filename, "DIC2DpairResults");

disp('File Saved - Finished')

