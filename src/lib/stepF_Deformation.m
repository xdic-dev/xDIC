function [filename] = stepF_Deformation(varargin)
    %% Parse input arguments
    p = inputParser;
    p.addParameter('basePath', []);
    p.addParameter('target_file', []);

    p.parse(varargin{:});
    basePath = p.Results.basePath;
    target_file = p.Results.target_file;

    %% Deformation reconstruction
    disp('Deformation reconstruction launch');
    tic;
    filename = step4_dic_rewrited(basePath, "cum", target_file);
    disp('step 4 : cum strains done');
    disp('--> STEP : 3D Deformation reconstruction completed');
    toc;
