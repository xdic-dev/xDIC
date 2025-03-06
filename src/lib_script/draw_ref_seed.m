function draw_ref_seed(base_parameters,varargin)

% [roi,~] = parseargpair(varargin,'roi');
p = inputParser;
p.addRequired('base_parameters');
p.addParameter('roi',[]);

p.parse(base_parameters,varargin{:});
roi = p.Results.roi;
[cam_ref_raw,~] = import_vid(base_parameters.baseDataPath,...
                    'subject',base_parameters.subject,...
                    'material',base_parameters.material,...
                    'trial',base_parameters.reftrial,...
                    'stereopair',base_parameters.stereopair,...
                    'phase',base_parameters.phase,...
                    'framejump',base_parameters.jump);
%%
message = {'Seed position not found!','Set the seed position in the reference trial and close when done'};
fig = figure; 
imshow(im_highlight(cam_ref_raw(:,:,1),roi,'factor',4),[])
seed_point = drawpoint(); 
seed_point = round(seed_point.Position);
title(message); 

%%
seed_file = base_parameters.seedfile; 
save(seed_file, "seed_point");

close(fig); 

end