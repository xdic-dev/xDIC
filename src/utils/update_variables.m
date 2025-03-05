if exist('subject', 'var') % if subject is defined subject_id = subject; end
    subject_id = subject;
end

if exist('phase', 'var') % if phase is defined phase_id = phase; end
    phase_id = phase;
end

if exist('material', 'var') % if material is defined material_id = material; end
    material_id = material;
end

if exist('nfcond', 'var') % if nfcond is defined nfcond_set = nfcond; end
    nfcond_set = nfcond;
end

if exist('spddxlcond', 'var') % if spddxlcond is defined spddxlcond_set = spddxlcond; end
    spddxlcond_set = spddxlcond;
end

if exist('calib_folder', 'var') % if calib_folder is defined calib_folder_set = calib_folder; end
    calib_folder_set = calib_folder;
end

if exist('ref_trial', 'var') % if ref_trial is defined ref_trial_id = ref_trial; end
    ref_trial_id = ref_trial;
end

if exist('frame_start', 'var') % if frame_start is defined idx_frame_start = idx_frame_start; end
    idx_frame_start = frame_start;
end

if exist('frame_end', 'var') % if frame_end is defined idx_frame_end = idx_frame_end; end
    idx_frame_end = frame_end;
end

if exist('jump', 'var') % if jump is defined frame_jump = jump; end
    frame_jump = jump;
end

if exist('visu', 'var') % if visu is defined showvisu = showvisu; end
    showvisu = visu;
end

if exist('debug', 'var') % if debug is defined debug_mode = debug; end
    debug_mode = debug;
end

material = frictional_conditions(material_id); % Material identifier
