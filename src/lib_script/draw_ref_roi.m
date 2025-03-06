function draw_ref_roi(base_parameters) 

%% read images of the reference trial
[cam_first_ref,~] = import_vid(base_parameters.baseDataPath,...
    'subject',base_parameters.subject,...
    'material',base_parameters.material,...
    'trial',base_parameters.reftrial,...
    'stereopair',base_parameters.stereopair,...
    'phase',base_parameters.phase,...
    'framejump',base_parameters.jump);
%% saturation of image
cam_first_ref = satur(cam_first_ref,'level',base_parameters.limit_grayscale);

%% start ncorr to create ROI
handles = ncorr;
handles.set_ref(cam_first_ref(:,:,1));
handles.set_cur(cam_first_ref(:,:,1));
disp('Draw ROI'); 

fig = uifigure;
message = {'ROI not found!','Draw ROI using ncorr and close when done'};
uialert(fig,message, ...
'Program Information','Icon','info','CloseFcn',@(h,e) close(fig)); 
uiwait(fig); 

refmask = handles.reference.roi.mask;
save(base_parameters.roifile, "refmask");

close(handles.handles_gui.figure);
disp('--> ROI cam first drawn'); 

end