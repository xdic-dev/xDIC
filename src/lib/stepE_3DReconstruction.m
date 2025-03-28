function filename = stepE_3DReconstruction(varargin)
    %% Parse input arguments
    p = inputParser;
    p.addParameter('basePath', []);
    p.addParameter('calibPath', []);
    p.addParameter('showvisu', 0);
    p.addParameter('title_ax', []);
    p.addParameter('pairOrder', [1 2]);
    p.addParameter('pairForced', false);

    p.parse(varargin{:});
    basePath = p.Results.basePath;
    calibPath = p.Results.calibPath;
    showvisu = p.Results.showvisu;
    title_ax = p.Results.title_ax;
    pairOrder = p.Results.pairOrder;
    pairForced = p.Results.pairForced;

    %% 3D reconstruction
    tic;
    fprintf("Starting 3D reconstruction:: stepE\n")
    filename = step3_dic_rewrited('basePath', basePath, ...
        'calibPath', calibPath, ...
        'showvisu', showvisu, ...
        'title_ax', title_ax, ...
        'pairOrder', pairOrder, ...
        'pairForced', pairForced);
    toc;
    disp('--> STEP : 3D reconstruction completed');
end
