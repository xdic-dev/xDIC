%% STEP 3: 3D reconstruction
% 3D reconstruction of points correlated with Ncorr using DLT 
% for running this step, you need to have DIC2DpairResults structures of
% the pairs you want to reconstruct from STEP2, as well as the DLT
% parameters of the cameras from STEP1.
function savefileName = step3_dic_rewrited(varargin)

%% Parse input arguments
p = inputParser;
p.addParameter('basePath',[]);
p.addParameter('calibPath',[]);
p.addParameter('showvisu',0);
p.addParameter('title_ax',[]);
p.addParameter('pairOrder',[]);
p.addParameter('pairForced',false);

p.parse(varargin{:});
basePath = p.Results.basePath;
calibPath = p.Results.calibPath;
showvisu = p.Results.showvisu;
title_ax = p.Results.title_ax;
pairOrder = p.Results.pairOrder;
pairForced = p.Results.pairForced;

fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

fprintf("Starting 3D reconstruction:: stepE, step3\n");

%% CHOOSE PATHS OPTIONS

% select DIC2DpairResults structures

structPaths{1} = fullfile(basePath,sprintf('myDIC2DpairResults_*C_1*.mat'));
structPaths{2} = fullfile(basePath,sprintf('myDIC2DpairResults_*C_4*.mat'));
nPairs=length(structPaths);
stitchButton = nPairs>1; 

fprintf("Starting 3D reconstruction:: stepE, step3,:: 1 \n");

indPairs=zeros(nPairs,2);
DIC2D=cell(nPairs,1);
for ip=1:nPairs
    S = dir(structPaths{ip});
    if isempty(S)
        error(sprintf('Error: pair#%d file not found.',ip));
    end
    DIC2D{ip}=load(fullfile(S.folder, S.name));
    DIC2D{ip}=DIC2D{ip}.DIC2DpairResults;
    indPairs(ip,:)=[DIC2D{ip}.nCamRef DIC2D{ip}.nCamDef];
end
indCams=unique(indPairs(:));
nCams=numel(indCams);

fprintf("Starting 3D reconstruction:: stepE, step3,:: 2 \n");

% select the folder where the DLT parameter files are stored
DLTparameterFolder = calibPath;
DLTstructAllCams=cell(nCams,1);
DLTpath=cell(nCams,1);
fprintf("Starting 3D reconstruction:: stepE, step3,:: 2, %s \n", DLTparameterFolder);
for ic=1:nCams
    fprintf("Starting 3D reconstruction:: stepE, step3,:: 2, %d \n", ic);
    Spath = fullfile(DLTparameterFolder,sprintf("*cam_%d.mat",indCams(ic)));
    fprintf("Starting 3D reconstruction:: stepE, step3,:: 2, %d, %s \n", ic, Spath );
    S = dir(Spath);
    DLTpath{ic} = fullfile(S.folder,S.name);
    DLTstructTemp=load(DLTpath{ic});
    DLTstructAllCams{ic}=DLTstructTemp.DLTstructCam;
end

fprintf("Starting 3D reconstruction:: stepE, step3,:: 3 \n");

% % remove distortion? If yes, choose the folder where the parameters are saved
% distortionRemovalButton = questdlg('Remove distortion?', 'Remove distortion?', 'Yes', 'No', 'Yes');
% switch distortionRemovalButton
%     case 'Yes'
%         distortionRemovalLogic=true(1);
%         distortionParametersPath=uigetdir([],'Select a folder where the cameraCBparameters are located');
%     case 'No'
%         distortionRemovalLogic=false(1);
% end 

% % save 3D-DIC results? choose save path and overwrite options
% [save3DDIClogic,savePath]=Qsave3DDICresults(structPaths);
save3DDIClogic = 1;
savePath = basePath; 
distortionRemovalLogic = 0; 
% %% if distortion removal is required, calculate the new DLT parameters
% switch distortionRemovalLogic
%     case 0
%         
%     case 1
%         distortionPar=cell(nCams,1);
%         distortionParPath=cell(nCams,1);
%         for ic=1:nCams
%             icam=indCams(ic);
%             P2Dtemp=DLTstructAllCams{ic}.imageCentroids;
%             columns=DLTstructAllCams{ic}.columns;
%             C3Dtrue=DLTstructAllCams{ic}.C3Dtrue;
%             % undistort image points
%             distortionParPath{ic}=fullfile(distortionParametersPath, ['cameraCBparameters_cam_' num2str(icam)]);
%             distortionParTemp=load(distortionParPath{ic});
%             distortionPar{ic}=distortionParTemp.cameraCBparameters.cameraParameters;
%             [P2Dtemp] = undistortPoints(P2Dtemp,distortionPar{ic});
%             % calculate new DLT parameters
%             C3Dtemp=reshape(C3Dtrue(:,columns,:),size(C3Dtrue,1)*numel(columns),3);
%             L=DLT11Calibration(P2Dtemp,C3Dtemp);
%             % assign to DLTstructPairs
%             DLTstructAllCams{ic}.DLTparams=L;
%         end
%     
% end

%% 3D reconstruction using Direct Linear Transformation
DIC3DAllPairsResults=cell(nPairs,1);

fprintf("Starting 3D reconstruction:: stepE, step3,:: 4 \n");

for ip=1:nPairs % loop over stereo pairs
%     hw = waitbar(0,['Reconstructing 3D points for pair ' num2str(ip) '...']);
    
    % create 3D-DIc results struct
    DIC3DpairResults=struct;
    
    % camera indices of current pair
    nCamRef=indPairs(ip,1);
    nCamDef=indPairs(ip,2);
    
    iCamRef=find(indCams==nCamRef);
    iCamDef=find(indCams==nCamDef);
    
    DIC3DpairResults.cameraPairInd=[nCamRef nCamDef];
    
    % DLT parameters
    L1=DLTstructAllCams{iCamRef}.DLTparams; % DLT parameters for 1st camera
    L2=DLTstructAllCams{iCamDef}.DLTparams; % DLT parameters for 2nd camera
    DIC3DpairResults.calibration.DLTpath{1}=DLTpath{iCamRef};
    DIC3DpairResults.calibration.DLTpath{2}=DLTpath{iCamDef};
    DIC3DpairResults.calibration.DLTparameters{1}=L1;
    DIC3DpairResults.calibration.DLTparameters{2}=L2;
            
    % extract information from 2d-dic results
    nImages=DIC2D{ip}.nImages;
    CorCoeff=DIC2D{ip}.CorCoeffVec;
    F=DIC2D{ip}.Faces;
    FC=DIC2D{ip}.FaceColors;

    DIC3DpairResults.Faces=F;
    DIC3DpairResults.FaceColors=FC;
    
    switch distortionRemovalLogic
        case 0 % if no distortion correction, use original points and set distortion model to 'none'
            Points=DIC2D{ip}.Points;
            DIC3DpairResults.distortionModel{1}='none';
            DIC3DpairResults.distortionModel{2}='none';
            DIC3DpairResults.distortionPath{1}='none';
            DIC3DpairResults.distortionPath{2}='none';
%         case 1 % if yes distortion correction, correct points and set distortion model
%             Points=DIC2D{ip}.Points;
%             DIC3DpairResults.distortionModel{1}=distortionPar{iCamRef};
%             DIC3DpairResults.distortionModel{2}=distortionPar{iCamDef};
%             DIC3DpairResults.distortionPath{1}=distortionParPath{iCamRef};
%             DIC3DpairResults.distortionPath{2}=distortionParPath{iCamDef};
%             step=min([40,size(Points{1},1)]);
%             % remove distortion from points
%             % first camera
%             for ii=1:nImages 
%                 %waitbar(ii/(6*nImages));
%                 
%                 % remove nan points
%                 PointsNaNlogic=any(isnan(Points{ii}),2);
%                 PointNoNan=Points{ii}(~PointsNaNlogic,:);
%                 distortionParNow=distortionPar{iCamRef};
%                 PointsTemp=zeros(size(PointNoNan));
%                 times=1:step:size(PointNoNan,1);
%                 % correct points in small sets because large sets are slow
%                 for jj=times(1:end-1)
%                     PointsTemp(jj:jj+step,:) = undistortPoints(PointNoNan(jj:jj+step,:),distortionParNow);
%                 end
%                 PointsTemp(times(end):end,:) = undistortPoints(PointNoNan(times(end):end,:),distortionParNow);
%                 Points{ii}(~PointsNaNlogic,:) = PointsTemp;                
%             end
%             % repeat for second camera
%             for ii=nImages+1:2*nImages 
%                 %waitbar(ii/(6*nImages));
%                 
%                 PointsNaNlogic=any(isnan(Points{ii}),2);
%                 PointNoNan=Points{ii}(~PointsNaNlogic,:);
%                 distortionParNow=distortionPar{iCamDef};
%                 PointsTemp=zeros(size(PointNoNan));
%                 times=1:step:size(PointNoNan,1);
%                 for jj=times(1:end-1)
%                     PointsTemp(jj:jj+step,:) = undistortPoints(PointNoNan(jj:jj+step,:),distortionParNow);
%                 end
%                 PointsTemp(times(end):end,:) = undistortPoints(PointNoNan(times(end):end,:),distortionParNow);
%                 Points{ii}(~PointsNaNlogic,:) = PointsTemp;
%             end 
            
    end

    % pre-allocate 3D-DIC result variables
    DIC3DpairResults.Points3D=cell(nImages,1);
    DIC3DpairResults.Disp.DispVec=cell(nImages,1);
    DIC3DpairResults.Disp.DispMgn=cell(nImages,1);
    DIC3DpairResults.FaceCentroids=cell(nImages,1);
    DIC3DpairResults.corrComb=cell(nImages,1);
    DIC3DpairResults.FaceCorrComb=cell(nImages,1);
    
    for ii=1:nImages % loop over images (time frames)
        %waitbar(1/3+ii/(3*nImages));
        % correlated points from 2 cameras
        P1=Points{ii};
        P2=Points{ii+nImages};
        
        % Solve the DLT system
        P3D=DLT11Reconstruction(P1,P2,L1,L2);
        DIC3DpairResults.Points3D{ii}=P3D;
        
        % Combined (worst) correlation coefficients
        DIC3DpairResults.corrComb{ii}=max([CorCoeff{ii} CorCoeff{ii+nImages}],[],2);
        % Face correlation coefficient (worst)
        DIC3DpairResults.FaceCorrComb{ii}=max(DIC3DpairResults.corrComb{ii}(F),[],2);
        
        % compute face centroids
        for iface=1:size(F,1)
            DIC3DpairResults.FaceCentroids{ii}(iface,:)=mean(P3D(F(iface,:),:));
        end

        % Compute displacements between frames (per point)
        DispVec=DIC3DpairResults.Points3D{ii}-DIC3DpairResults.Points3D{1};
        DIC3DpairResults.Disp.DispVec{ii}=DispVec;
        DIC3DpairResults.Disp.DispMgn{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
        
    end
    
    % put all pairs in a cell array
    DIC3DAllPairsResults{ip}=DIC3DpairResults;
    
%      delete(hw) 
end


%% Stitch pairs
% if nPairs>1
%     % Stitch pairs? If yes, selesct which pairs to stitch and in which order
%     stitchButton = questdlg('Stitch surfaces together?', 'Stitch surfaces together?', 'Yes', 'No', 'Yes');
%     switch stitchButton
%         case 'Yes'
%             stitchButton=true(1);
%         case 'No'
%             stitchButton=false(1);
%     end
% else
%     stitchButton=false(1);
% end
if stitchButton
    % stitch or stitch+append
    %anim8_DIC3D_reconstructedPairs_faceMeasure(DIC3DAllPairsResults,'pairInd');
%     answer = inputdlg({sprintf('Enter the indices of camera-pairs to stitch to each other (in the order they should be stitched):\nSurfaces deleted from the list will not be stitched, but will be included in the results')},'Input',[1,100],{mat2str(1:nPairs)});
    pairIndList=pairOrder;
    [DIC3Dcombined]= DIC3DsurfaceStitch(DIC3DAllPairsResults,pairIndList,...
        'showvisu',showvisu,'title_ax',title_ax,'pairForced',pairForced);
else
    % only append
    [DIC3Dcombined]= DIC3DsurfaceStitch(DIC3DAllPairsResults,[]);
end

% add all information to DIC3Dcombined structure
for ipair=1:nPairs
    % add distortion model information to DIC3Dcombined
    DIC3Dcombined.distortion.distortionPath{ipair,1}=DIC3DAllPairsResults{ipair}.distortionPath{1};
    DIC3Dcombined.distortion.distortionPath{ipair,2}=DIC3DAllPairsResults{ipair}.distortionPath{2};
    DIC3Dcombined.distortion.distortionModel{ipair,1}=DIC3DAllPairsResults{ipair}.distortionModel{1};
    DIC3Dcombined.distortion.distortionModel{ipair,2}=DIC3DAllPairsResults{ipair}.distortionModel{2};
    % add calibration information to DIC3Dcombined
    DIC3Dcombined.calibration.DLTpath{ipair,1}=DIC3DAllPairsResults{ipair}.calibration.DLTpath{1};
    DIC3Dcombined.calibration.DLTpath{ipair,2}=DIC3DAllPairsResults{ipair}.calibration.DLTpath{2};
    DIC3Dcombined.calibration.DLTparameters{ipair,1}=DIC3DAllPairsResults{ipair}.calibration.DLTparameters{1};
    DIC3Dcombined.calibration.DLTparameters{ipair,2}=DIC3DAllPairsResults{ipair}.calibration.DLTparameters{2};
    % add 2D-DIC data to DIC3Dcombined
    DIC3Dcombined.DIC2Dinfo{ipair}=DIC2D{ipair};
    % add 3D-DIC data to DIC3Dcombined
    DIC3Dcombined.AllPairsResults=DIC3DAllPairsResults;
end

%% save 

if save3DDIClogic
    if stitchButton
        savefileName = ['DIC3Dcombined_' num2str(nPairs) 'Pairs_stitched.mat'];
        saveName=fullfile(savePath,savefileName);
    else
        savefileName = ['DIC3Dcombined_' num2str(nPairs) 'Pairs_stitched.mat'];
        saveName=fullfile(savePath,savefileName);
    end
    icount=1;
    while exist(saveName,'file')
        if stitchButton
            savefileName = ['DIC3Dcombined_' num2str(nPairs) 'Pairs_stitched_v' num2str(icount) '.mat'];
            saveName=fullfile(savePath,savefileName);
        else
            savefileName = ['DIC3Dcombined_' num2str(nPairs) 'Pairs_stitched_v' num2str(icount) '.mat'];
            saveName=fullfile(savePath,savefileName);
        end
        icount=icount+1;
    end
    save(saveName,'DIC3Dcombined','-v7.3');
end

% %% finish
% h=msgbox('STEP3 is completed');
% h.CurrentAxes.Children.FontSize=11;
% 
% set(0, 'DefaultUIControlFontSize', fs);

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>
end