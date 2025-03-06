function saveName = step4_dic_rewrited(basePath,deftype,reconstruct_file)
%% STEP 4: Post Processing
% Upload 3D reconstructed points (either as seperate pairs or a stitched surface)
% and calculate diplacements, deformations, and strains.

%%
% clearvars; close all
switch deftype
    case "rate"
        cum = 0;  
    case "cum"
        cum = 1;  
end
cum = 1; 
fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

%% CHOOSE PATHS OPTIONS

% select DIC3DpairResults structures
% PathInitial=pwd;
% [file,path] = uigetfile(PathInitial,'Select a 3D-DIC results structure (stitched or unstitched)');
path = basePath; 
if isempty(reconstruct_file) 
    file = "DIC3Dcombined_1Pairs_stitched.mat"; 
else 
    file = reconstruct_file; 
end
DIC3D=load(fullfile(path,file));
DIC3Dname=fieldnames(DIC3D);
DIC3D=DIC3D.(DIC3Dname{1});

% save 3D-DIC post processing results? choose save path and overwrite options
% [save3DDIClogic,savePath]=Qsave3DDICPPresults(path);
save3DDIClogic = 1;
savePath = basePath;

%% 3D reconstruction using Direct Linear Transformation

nImages= numel(DIC3D.Points3D);

% pre-allocate 3D-DIC result variables
% DIC3D.Points3D_ARBM=cell(1,nImages);
DIC3D.Disp.DispVec=cell(1,nImages);
DIC3D.Disp.DispMgn=cell(1,nImages);
% DIC3D.Disp.DispVec_ARBM=cell(1,nImages);
% DIC3D.Disp.DispMgn_ARBM=cell(1,nImages);
DIC3D.FaceCentroids=cell(1,nImages);
% DIC3D.FaceCentroids_ARBM=cell(1,nImages);
DIC3D.FaceCorrComb=cell(1,nImages);
DIC3D.FaceIsoInd=cell(1,nImages);
% DIC3D.RBM.RotMat=cell(1,nImages);
% DIC3D.RBM.TransVec=cell(1,nImages);

F=DIC3D.Faces;
% hw = waitbar(0,'Calculating displacements and rigid body motion');
for ii=1:nImages % loop over images (time frames)
%     waitbar(ii/(nImages));
    
    % Face correlation coefficient (worst)
    DIC3D.FaceCorrComb{ii}=max(DIC3D.corrComb{ii}(F),[],2);
    
    % compute face centroids
    for iface=1:size(F,1)
        DIC3D.FaceCentroids{ii}(iface,:)=mean(DIC3D.Points3D{ii}(F(iface,:),:));
    end
    
    % Compute displacements between frames (per point)
    DispVec=DIC3D.Points3D{ii}-DIC3D.Points3D{1};
    DIC3D.Disp.DispVec{ii}=DispVec;
    DIC3D.Disp.DispMgn{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
    
%     % Compute rigid body transformation between point clouds
%     [RotMat,TransVec,Points3D_ARBM]=rigidTransformation(DIC3D.Points3D{ii},DIC3D.Points3D{1});
%     DIC3D.RBM.RotMat{ii}=RotMat;
%     DIC3D.RBM.TransVec{ii}=TransVec;
%     DIC3D.Points3D_ARBM{ii}=Points3D_ARBM;
%     
%     % Compute displacements between sets - after RBM
%     DispVec=DIC3D.Points3D_ARBM{ii}-DIC3D.Points3D_ARBM{1};
%     DIC3D.Disp.DispVec_ARBM{ii}=DispVec;
%     DIC3D.Disp.DispMgn_ARBM{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
%     
%     % compute face centroids - after transformation
%     for iface=1:size(F,1)
%         DIC3D.FaceCentroids_ARBM{ii}(iface,:)=mean(Points3D_ARBM(F(iface,:),:));
%     end
    
end
% delete(hw);

% %filtering process before deformation computation
% freqFilt = 3; 
% actual_FS_vid = 50/5; 
% for kk = 1:3 
%     DIC3D.Disp.DispVec = filter3Ddeform_time(DIC3D.Disp.DispVec,...
%         'freqFilt',freqFilt,...
%         'freqAcq',actual_FS_vid,...
%         'deformCol',kk);
% end
% DIC3D.Disp.DispMgn = filter3Ddeform_time(DIC3D.Disp.DispMgn,...
%     'freqFilt',freqFilt,...
%     'freqAcq',actual_FS_vid);

% compute deformation and strains (per triangular face)
deformationStruct=triSurfaceDeformation_rewrited(F,DIC3D.Points3D{1},DIC3D.Points3D,cum);
DIC3D.Deform=deformationStruct;

% deformationStruct_ARBM=triSurfaceDeformation(F,DIC3D.Points3D_ARBM{1},DIC3D.Points3D_ARBM);
% DIC3D.Deform_ARBM=deformationStruct_ARBM;

% compute triangle regularity (isotropy index)
for ii=1:nImages
    [FisoInd]=faceIsotropyIndex(F,DIC3D.Points3D{ii});
    DIC3D.FaceIsoInd{ii}=FisoInd;
end

DIC3DPPresults=DIC3D;


%% save results
nPairs=size(DIC3DPPresults.pairIndices,1);
if save3DDIClogic
    saveName=fullfile(savePath, "DIC3DPPresults_"+num2str(nPairs)+"Pairs_"+deftype+"_v1.mat");
    icount=1;
    while exist(saveName,'file')
        saveName=fullfile(savePath, "DIC3DPPresults_"+num2str(nPairs)+"Pairs_"+deftype+"_v"+num2str(icount+1)+".mat");
        icount=icount+1;
    end
    save(saveName,'DIC3DPPresults','-v7.3');
end
set(0, 'DefaultUIControlFontSize', fs);

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