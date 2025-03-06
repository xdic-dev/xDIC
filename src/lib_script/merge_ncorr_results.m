function ncorr_modif = merge_ncorr_results(ncorr1,ncorr2)
%%% Insert at the beginning ncorr1 and push the second ncorr analysis. 
    
    %Ncorr length
    N1 = length(ncorr1.current_save);
    N2 = length(ncorr2.current_save);
    
    %Initialisation of new ncorr
    ncorr_modif = ncorr2;
    
    %field names of structures
    disp_field = fieldnames(ncorr_modif.data_dic_save.displacements);
    ref_field = fieldnames(ncorr_modif.reference_save);
    cur_field = fieldnames(ncorr_modif.current_save);
    
    %tail - ncorr2
    for ii = 2:N2
        pp = N2+2-ii; 
        for jj = 1:length(disp_field)
            if (~strcmp(disp_field{jj},'roi_dic')&&~strcmp(disp_field{jj},'roi_ref_formatted')&&~strcmp(disp_field{jj},'roi_cur_formatted'))
                ncorr_modif.data_dic_save.displacements((N1-1)+pp).(disp_field{jj}) = ...
                    ncorr1.data_dic_save.displacements(N1).(disp_field{jj})+...
                    ncorr2.data_dic_save.displacements(pp).(disp_field{jj});
            else
                if pp > N2-N1+1
                    %mask
                    ncorr_modif.data_dic_save.displacements((N1-1)+pp).(disp_field{jj}).mask = ...
                        ncorr2.data_dic_save.displacements(pp).(disp_field{jj}).mask;
                end
            end
        end
        for kk = 1:length(cur_field)
            ncorr_modif.current_save((N1-1)+pp).(cur_field{kk}) = ...
                ncorr2.current_save(pp).(cur_field{kk});
        end
    end
    %head - ncorr1
    for ii = 1:N1
        pp = ii; 
        for jj = 1:length(disp_field)
            if ~strcmp(disp_field{jj},'roi_dic')
                ncorr_modif.data_dic_save.displacements(pp).(disp_field{jj}) = ...
                    ncorr1.data_dic_save.displacements(pp).(disp_field{jj});
            else
%                 %mask
%                 ncorr_modif.data_dic_save.displacements(ii).(disp_field{jj}).mask = ...
%                     ncorr1.data_dic_save.displacements(ii).(disp_field{jj}).mask;
            end
        end
        for kk = 1:length(cur_field)
            ncorr_modif.current_save(pp).(cur_field{kk}) = ...
                ncorr1.current_save(pp).(cur_field{kk});
        end
    end
    
    %reference - ncorr1
    for jj = 1:length(ref_field)
        ncorr_modif.reference_save.(ref_field{jj}) = ncorr1.reference_save.(ref_field{jj});
    end
        
end