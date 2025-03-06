
function data=remove_highcorr_level(data,level)
fn = fieldnames(data.Deform); 
Ndeform = numel(fn); 
Nframe = length(data.corrComb);
mask_highcorr = cell(Nframe,1); 
for jj = 1:Ndeform
    deformjj = fn{jj}; 
    for ii = 1:Nframe
        if jj == 1
            mask_highcorr{ii} = data.corrComb{ii}>level;
        end
        data.(deformjj){ii}(mask_highcorr{ii}) = NaN; 
    end
end
