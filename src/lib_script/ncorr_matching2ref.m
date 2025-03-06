function ncorr_matching2ref(im2match,base_parameters,ncorr_parameters)

disp('Ncorr matching run'); 

%% Load images
[cam_first_raw,~] = import_vid(base_parameters.baseDataPath,...
                    'subject',base_parameters.subject,...
                    'material',base_parameters.material,...
                    'trial',base_parameters.reftrial,...
                    'stereopair',base_parameters.stereopair,...
                    'phase',base_parameters.phase,...
                    'idxstart_set',base_parameters.idxstart_set,...
                    'idxend_set',base_parameters.idxstart_set+1,...
                    'framejump',base_parameters.jump);
%
cam_first = cam_first_raw;

%% REFERENCE ROI 
refmask = load(base_parameters.roifile); 
refmask = refmask.refmask;

%% Run ncorr 
%Choose automatic or semi automatic analysis

if exist(base_parameters.seedfile,'file') 
    seed_point = load(base_parameters.seedfile); 
    seed_point = seed_point.seed_point;
    showGui = false; 
else 
    seed_point = []; 
    showGui = false; 
end

siz = size(im2match); 
input2 = zeros(siz(1),siz(2),2,'uint8'); 
input2(:,:,1) = im2match; 
input2(:,:,2) = cam_first(:,:,1); 

handles = ncorr_abr(showGui,seed_point);
handles.set_ref(cam_first(:,:,1));
handles.set_cur(Myarray2cell(input2));
handles.set_roi_ref(refmask);

if (true)
    handles.set_dic_parameters_manually( ...
        ncorr_parameters.type, ...
        ncorr_parameters.radius, ...
        ncorr_parameters.spacing, ...
        ncorr_parameters.cutoff_diffnorm, ...
        ncorr_parameters.cutoff_iteration,...
        ncorr_parameters.total_threads, ...
        ncorr_parameters.stepanalysis_params)

    % refresh the handle
    handles.abr_refresh()
    
    % trigger DIC
    handles.callback_topmenu_dic();

    % trigger Displacement
    handles.callback_topmenu_formatdisp();

end

%% Save ncorr results
file = base_parameters.matchingfile; 
save_dicdata(handles, file); 
close(handles.handles_gui.figure);
disp('--> Ncorr matching done');
disp('--> Ncorr analysis completed');


end


% U_mapped = imresize(matching.data_dic_save.displacements(1).plot_u_ref_formatted,size(refmask),'nearest'); 
% V_mapped = imresize(matching.data_dic_save.displacements(1).plot_v_ref_formatted,size(refmask),'nearest'); 
% 
% [U_mapped_row, U_mapped_col, U_mapped_value] = find(U_mapped);
% [V_mapped_row, V_mapped_col, V_mapped_value] = find(V_mapped);
% 
% newmask = false(size(U_mapped)); 
% newmask(round(U_mapped_row+U_mapped_value), round(U_mapped_col+V_mapped_value)) = true; 
% % newmask(~refmask) = false; 
% oldmask = refmask; 
% refmask = newmask; 




