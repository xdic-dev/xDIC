function [varargout]=anim8_DIC3DPP_faceMeasure_rewrited(data,faceMeasureString,varargin)
%% function for plotting 3D-DIC results of face measures as color + direction as arrows in STEP4.
% plotting 3D surfaces from camera pairs, animation changing
% with time, and the faces colored according to faceMeasureString 
% this function is called in plotMultiDICPairResults
%
% Options:
% anim8_DIC3DPP_faceMeasureDirection(DIC_3Dallpairs_results,faceMeasureString)
% anim8_DIC3DPP_faceMeasureDirection(DIC_3Dallpairs_results,faceMeasureString,optStruct)
% 
% Inputs:
% * DIC3DAllPairsResults
% * faceMeasureString: can be any of the following:
%   'Epc1','Epc2','epc1','epc2','Lamda1','Lamda2', or a combination of two  in a cell array to plot side by side
% * optStruct: optional structure for plotting options which may include any of the following fields:
%   - smoothLogic: logical variable for smoothing (true)/not smoothing (false) the face measure 
%   - FaceAlpha: transparacy of the faces (scalar between 0 and 1, where zero is transparent and 1 is opaque) 
%   - colorBarLimits: a 2x1 scalar vector for the colobar limits. if not set, it's automatic
%   - dataLimits: a 2x1 scalar vector for the data limits of the face measure. if a face measure is outside these limits, it is set to NaN. if not set no face is set to NaN
%   - colorMap
%   - zDirection: 1 for z up and -1 for z down
%   - lineColor: line color for the mesh. can be for example 'b','k','none',etc...
%   - supTitleString=faceMeasureString;

%% Assign plot options
Narg=numel(varargin);
switch Narg
    case 1
        optStruct=varargin{1};
    case 0
        optStruct=struct;
    otherwise
        ('wrong number of input arguments');
end

% complete the struct fields
if ~isfield(optStruct,'smoothSpaceLogic')
    optStruct.smoothSpaceLogic=0;
    fprintf('no smoothSpaceLogic\n'); 
end
if ~isfield(optStruct,'FaceAlpha')
    optStruct.FaceAlpha=1;
end
if ~isfield(optStruct,'dataLimits')
    optStruct.dataLimits=[-inf inf];
end
if ~isfield(optStruct,'zDirection') % 1 or -1
    optStruct.zDirection=1;
end
if ~isfield(optStruct,'lineColor') % 'none' or 'k'
    optStruct.lineColor='none';
end
if ~isfield(optStruct,'maxCorrCoeff')
    optStruct.maxCorrCoeff=[];
end
if ~isfield(optStruct,'quiverScaleFactor')
    optStruct.quiverScaleFactor=20;
end
if ~isfield(optStruct,'maxabs_tot')
    optStruct.maxabs_tot=[];
end
if ~isfield(optStruct,'maxabs_rate')
    optStruct.maxabs_rate=[];
end
if ~isfield(optStruct,'showvisu')
    optStruct.showvisu=1;
end
if ~isfield(optStruct,'deftype')
    optStruct.deftype='cum';
end
if ~isfield(optStruct,'gapLogic')
    optStruct.gapLogic=0;
end
if ~isfield(optStruct,'showRobotLogic')
    optStruct.showRobotLogic=0;
end
if ~isfield(optStruct,'showImgLogic')
    optStruct.showImgLogic=0;
end
if ~isfield(optStruct,'contactAreaLogic')
    optStruct.contactAreaLogic = 0; 
end

%%
nFrames=numel(data.Points3D);
[xl,yl,zl]=axesLimits(data.Points3D);
meanEdgeLength=mean(patchEdgeLengths(data.Faces,data.Points3D{1}),'omitnan');
switch optStruct.deftype
    case 'cum'
        nDefType = 1; 
    case 'rate'
        nDefType = 1; 
    case 'both'
        nDefType = 2; 
end
if optStruct.showRobotLogic
    data_robot = optStruct.data_robot;
end
if optStruct.showImgLogic 
    data_img = optStruct.data_img; 
end

%% Assign the correct face measure into FC
if iscell(faceMeasureString)
    faceMeasureCell=faceMeasureString;
    nStrains=numel(faceMeasureString);
    FC=cell(nFrames,nStrains,nDefType);
    D=cell(nFrames,nStrains,nDefType);
    Ds=cell(nFrames,nStrains,nDefType);
    Vc=cell(nFrames,nStrains,nDefType);
    
elseif ischar(faceMeasureString)
    nStrains=1;
    faceMeasureCell=cell(1);
    faceMeasureCell{1}=faceMeasureString;
    FC=cell(nFrames,nDefType);
    D=cell(nFrames,nDefType);
    Ds=cell(nFrames,nDefType);
    Vc=cell(nFrames,nDefType);
    
else
    error('wrong face measure (second input variable)');
end

nvecplot = 0; 
directionStringCell = cell(nStrains,1);
colorBarLimitsCell = cell(nStrains,1);
% threshvalCell = cell(nStrains,1);
F=data.Faces;
for is=1:nStrains
    faceMeasureCell_is = faceMeasureCell{is}; 
    if strcmp(faceMeasureCell_is,'DispX')
        dispNow=data.Disp.DispVec{1}(:,1); % point measure
        z = zeros(nFrames,length(mean(dispNow(F),2)));
        FCmat = z;
        for it = 1:nFrames
            dispNow=data.Disp.DispVec{it}(:,1); % point measure
            FCmat(it,:)=mean(dispNow(F),2);
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    elseif strcmp(faceMeasureCell_is,'DispY')
        dispNow=data.Disp.DispVec{1}(:,2); % point measure
        z = zeros(nFrames,length(mean(dispNow(F),2)));
        FCmat = z;
        for it = 1:nFrames
            dispNow=data.Disp.DispVec{it}(:,2); % point measure
            FCmat(it,:)=mean(dispNow(F),2);
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    elseif strcmp(faceMeasureCell_is,'DispXY')
        dispNow=data.Disp.DispVec{1}(:,1); % point measure
        z = zeros(nFrames,length(mean(dispNow(F),2)));
        FCmat = z;
        for it = 1:nFrames
            dispNow=data.Disp.DispVec{it}(:,1); % point measure
            dispNow=dispNow+data.Disp.DispVec{it}(:,2); % point measure
            FCmat(it,:)=mean(dispNow(F),2);
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    elseif strcmp(faceMeasureCell_is,'DispZ')
        dispNow=data.Disp.DispVec{1}(:,3); % point measure
        z = zeros(nFrames,length(mean(dispNow(F),2)));
        FCmat = z;
        for it = 1:nFrames
            dispNow=data.Disp.DispVec{it}(:,3); % point measure
            FCmat(it,:)=mean(dispNow(F),2);
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;%-FCmat(24,:);
    elseif strcmp(faceMeasureCell_is,'DispMgn')
        dispNow=data.Disp.DispMgn{1}; % point measure
        z = zeros(nFrames,length(mean(dispNow(F),2)));
        FCmat = z;
        for it = 1:nFrames
            dispNow=data.Disp.DispMgn{it}; % point measure
            FCmat(it,:)=mean(dispNow(F),2);
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    elseif strcmp(faceMeasureCell_is,'corrComb')
        z = zeros(nFrames,length(data.(faceMeasureCell_is){1}));
        FCmat = z;
        for it = 1:nFrames
            FCmat(it,:)=data.(faceMeasureCell_is){it};
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    elseif isfield(data,faceMeasureCell_is) 
        z = zeros(nFrames,length(data.(faceMeasureCell_is)));
        FCmat = z;
        for it = 1:nFrames
            FCmat(it,:)=data.(faceMeasureCell_is);
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    else
        z = zeros(nFrames,length(data.Deform.(faceMeasureCell_is){1}));
        FCmat = z;
        for it = 1:nFrames
            FCmat(it,:)=data.Deform.(faceMeasureCell_is){it};
        end
        FCmat_rate = [z(1,:);diff(FCmat,1,1)];
        FCmat_cum = FCmat;
    end
    switch faceMeasureCell_is
        case 'Epc1'
            directionStringCell{is}='Epc1vecCur';
        case 'Epc2'
            directionStringCell{is}='Epc2vecCur';
        case 'epc1'
            directionStringCell{is}='epc1vec';
        case 'epc2'
            directionStringCell{is}='epc2vec';
        case 'Lamda1'
            directionStringCell{is}='Epc1vecCur';
        case 'Lamda2'
            directionStringCell{is}='Epc2vecCur';
        otherwise
            directionStringCell{is}=[];
    end
    maxval = zeros(nDefType,1);
    threshval = zeros(nFrames,nStrains,nDefType);
    if ~isempty(directionStringCell{is})
        nvecplot = nvecplot+1;
        for it=1:nFrames
            switch optStruct.deftype
                case 'cum'
                    FC{it,is,1} = FCmat_cum(it,:)';
                case 'rate'
                    FC{it,is,1} = FCmat_rate(it,:)';
                case 'both'
                    FC{it,is,2} = FCmat_cum(it,:)';
                    FC{it,is,1} = FCmat_rate(it,:)';
            end
            % direction (unit vector)
            D{it,is}=data.Deform.(directionStringCell{is}){it};
            
            switch faceMeasureCell_is
                case {'Epc1','Epc2','epc1','epc2'}
                    % direction with magnitude (scaled vector)
                    Ds{it,is}=optStruct.quiverScaleFactor*FC{it,is,1}.*D{it,is};
                    DsLengths=sqrt((sum(Ds{it,is}.^2,2)));
                    LogicTooLong=DsLengths>meanEdgeLength;
                    Ds{it,is}(LogicTooLong,:)=meanEdgeLength*Ds{it,is}(LogicTooLong,:)./DsLengths(LogicTooLong);
                case {'Lamda1','Lamda2'}
                    % direction with magnitude (scaled vector)
                    Ds{it,is}=optStruct.quiverScaleFactor*(FC{it,is,1}-1).*D{it,is};
                    DsLengths=sqrt((sum(Ds{it,is}.^2,2)));
                    LogicTooLong=DsLengths>meanEdgeLength;
                    Ds{it,is}(LogicTooLong,:)=meanEdgeLength*Ds{it,is}(LogicTooLong,:)./DsLengths(LogicTooLong);
            end
            Vc{it,is}=data.FaceCentroids{it}-.5*Ds{it,is};
            
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=data.FaceCorrComb{it};
                for idef = 1:nDefType
                    FC{it,is,idef}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
                end
                D{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
                Vc{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
            if ~isempty(optStruct.maxabs_tot)
                mask = abs(FC{it,is})>optStruct.maxabs_tot;
                for idef = 1:nDefType
                    FC{it,is,idef}(mask,:)=NaN;
                end
            end
            for idef = 1:nDefType
                maxval_iter = prctile(abs(FC{it,is,idef}),95);
                if maxval_iter>maxval(idef)
                    maxval(idef) = maxval_iter;
                end
                threshval(it,is,idef) = prctile(abs(FC{it,is,idef}),99);
            end
        end
    else %without vector
        for it=1:nFrames
            switch optStruct.deftype
                case 'cum'
                    FC{it,is,1} = FCmat_cum(it,:)';
                case 'rate'
                    FC{it,is,1} = FCmat_rate(it,:)';
                case 'both'
                    FC{it,is,2} = FCmat_cum(it,:)';
                    FC{it,is,1} = FCmat_rate(it,:)';
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=data.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
            if ~isempty(optStruct.maxabs_tot)
                mask = abs(FC{it,is})>optStruct.maxabs_tot;
                FC{it,is}(mask,:)=NaN;
            end
            for idef = 1:nDefType
                maxval_iter = prctile(abs(FC{it,is,idef}),95);
                if maxval_iter>maxval(idef)
                    maxval(idef) = maxval_iter;
                end
                threshval(it,is,idef) = prctile(abs(FC{it,is,idef}),99.5);
            end
        end
    end
    if ~isfield(optStruct,'colorBarLimits')
        for idef = 1:nDefType
            colorBarLimitsCell{is,idef}=[-round(maxval(idef)*1e3)/1e3 round(maxval(idef)*1e3)/1e3];
        end
    end
end

%% PLOT INITIALISATION 
% initialisation of figure properties
animStruct=struct;
hf=cFigure;
set(hf,'visible',optStruct.showvisu); 
hf.Units='normalized';
if isfield(optStruct,'FigPosition')
    hf.OuterPosition=optStruct.FigPosition;
else
    hf.OuterPosition=[.05 .05 .9 .9];
end
hf.Units='pixels';

% Issue when plotting only one strain type with forces plot
%--------------------------------------------------------------------------
% Initialisation of subplot handles
%--------------------------------------------------------------------------
Ncol = nStrains+optStruct.showRobotLogic+optStruct.showImgLogic; 
Nrow = nDefType; 
mergeElement = []; 
if showImgLogic
    mergeElement = [Ncol,Ncol*Nrow];
end
h_ax = subplot_ax(Nrow,Ncol,'merge',{mergeElement}); hold(h_ax,'on'); Nax = length(h_ax); 
hp = gobjects(Nax,1);
hp2 = gobjects(Nax,1);
hb = gobjects(Nax,1);
hq = gobjects(Nax,1);

it=1;
Fnow=data.Faces; %Face index 
for is=1:nStrains
    faceMeasureCell_is = faceMeasureCell{is}; 
    for idef = 1:nDefType
        ax_idx = is+(idef-1)*Ncol; 
        ax_cur = h_ax(ax_idx);

        %Settings for current axes
        set(hf,'CurrentAxes',ax_cur); 
        axisGeom(ax_cur);
        ax_cur.CameraUpVector=[0 0 optStruct.zDirection];
        colorbar(ax_cur);
        caxis(ax_cur,colorBarLimitsCell{is,idef});
%         title(ax_cur,faceMeasureCell_is);

        %current element
        FCnow=FC{it,is,idef}; %Value of each face
        Pnow=data.Points3D{it}; %Position of each vertex
        
        if optStruct.contactAreaLogic
            Fcentroidnow=data.FaceCentroids{it}; %Position of each face centroids
            Fcentroidlast=data.FaceCentroids{end}; %point cloud from last iteration
        end
        
        if ~isempty(directionStringCell{is})
            Vnow=Vc{it,is};
            Dnow=Ds{it,is};
        end
        
        %remove of unwanted elements
        FCnow(FCnow<optStruct.dataLimits(1))=NaN;
        FCnow(FCnow>optStruct.dataLimits(2))=NaN;
        FCnow(abs(FCnow)>threshval(it,is,idef))=NaN; %Maximum value removal
        if optStruct.gapLogic == 1
            %remove mesh between pairs : too noisy 
            FCnow(data.FacePairInds==3)=NaN; 
        end

        %spatial filtering 
        if optStruct.smoothSpaceLogic
            FCnow=patchSmoothFaceMeasure(Fnow,Pnow,FCnow,optStruct);
        end

        %contact area
        if optStruct.contactAreaLogic
            [CAnow,~,threshold_lastmap] = find_contact_area(Fcentroidnow,Fcentroidlast,'guess_threshold',-1); 
        end
        
        %plot creation 
        hp(ax_idx)=gpatch(Fnow,Pnow,FCnow,optStruct.lineColor,optStruct.FaceAlpha); hold on
        hp2(ax_idx)=gpatch(Fnow,Pnow.*0.92+[0 0 0.3]); hold on
        if optStruct.contactAreaLogic
            hb(ax_idx)=plot3(CAnow(:,1),CAnow(:,2),CAnow(:,3),'k-','LineWidth',2); hold on
        end
        if ~isempty(directionStringCell{is})
            hq(ax_idx)=quiver3(Vnow(:,1),Vnow(:,2),Vnow(:,3),Dnow(:,1),Dnow(:,2),Dnow(:,3),0,'Color',.2*[1 1 1],'ShowArrowHead','off','AutoScale','off'); hold on;
        end

        %limitation of axes
        ax_cur.XLim = xl; ax_cur.YLim = yl; ax_cur.ZLim = zl;
    end
end

% ROBOT RAW KINEMATICAL AND DYNAMICAL DATA  
if optStruct.showRobotLogic
   LINE_WIDTH = 3; 
   HRULE_WIDTH = 5; 
    %displacement 
    ax_idx = (nStrains+1)+Ncol*0; 
    set(hf,'CurrentAxes',h_ax(ax_idx)); 

    timestart = data_robot.time(optStruct.frame_idx(1));
    timeend = data_robot.time(optStruct.frame_idx(end));
    timeframe = optStruct.frame_idx(1):optStruct.frame_idx(end);
    dispX = data_robot.position_robot(timeframe,1)-data_robot.position_robot(optStruct.frame_idx(1),1);
    dispY = data_robot.position_robot(timeframe,2)-data_robot.position_robot(optStruct.frame_idx(1),2);
    disp = sqrt(dispX.^2+dispY.^2);
    maxvalrule(1) = ceil(max(disp)); 
    %
    plot(h_ax(ax_idx),data_robot.time(timeframe)-timestart,disp,'k-','LineWidth',LINE_WIDTH);
    %
    h_rule(1)=plot(h_ax(ax_idx),[0 0],[0 maxvalrule(1)],'r-','LineWidth',HRULE_WIDTH);
    %
    h_ax(ax_idx).XLabel.String = 'Time (s)'; h_ax(ax_idx).YLabel.String = 'plate disp (mm)';
    h_ax(ax_idx).XLim = [0 timeend-timestart]; h_ax(ax_idx).YLim = [0 maxvalrule(1)];

    %force 
    ax_idx = (nStrains+1)+Ncol*1; 
    set(hf,'CurrentAxes',h_ax(ax_idx)); 

    timestart = data_robot.time(optStruct.frame_idx(1));
    timeend = data_robot.time(optStruct.frame_idx(end));
    maxvalrule(2) = ceil(max(data_robot.tf(optStruct.frame_idx))); 
    if maxvalrule(2)<ceil(max(data_robot.nf(optStruct.frame_idx)))
        maxvalrule(2) = ceil(max(data_robot.nf(optStruct.frame_idx)));
    end
    timeframe = optStruct.frame_idx(1):optStruct.frame_idx(end);
    %
    plot(h_ax(ax_idx),data_robot.time(timeframe)-timestart,data_robot.nf(timeframe),'k-','LineWidth',LINE_WIDTH);
    plot(h_ax(ax_idx),data_robot.time(timeframe)-timestart,data_robot.tf(timeframe),'k-','LineWidth',LINE_WIDTH,'Color',[1 1 1].*0.5);
    %
    h_rule(2)=plot(h_ax(ax_idx),[0 0],[0 maxvalrule(2)],'r-','LineWidth',HRULE_WIDTH);
    %
    h_ax(ax_idx).XLabel.String = 'Time (s)'; h_ax(ax_idx).YLabel.String = 'Force (N)';
    h_ax(ax_idx).XLim = [0 timeend-timestart]; h_ax(ax_idx).YLim = [0 maxvalrule(2)];
    % 
    legend(h_ax(ax_idx),{'NF','TF'}); 
end

% RAW IMAGES FROM ONE CAMERA 
if optStruct.showImgLogic 
    ax_idx = (nStrains+2); 
    set(hf,'CurrentAxes',h_ax(ax_idx)); 
    himg = imshow(data_img(:,:,it),[]);
end

%% FILL PLOT DATA 
animStruct.Time=1:nFrames;
animStruct.Handles=cell(1,nFrames);
animStruct.Props=cell(1,nFrames);
animStruct.Set=cell(1,nFrames);

for it=1:nFrames
    
    nProps = 12; 
    animStruct.Handles{it}=[];
    animStruct.Props{it}=cell(1,nProps*nStrains*nDefType);
    animStruct.Set{it}=cell(1,nProps*nStrains*nDefType);
    
    ii_anim = 0; 
    for is=1:nStrains
        for idef = 1:nDefType
            ax_idx = is+(idef-1)*Ncol;
            ax_cur = h_ax(ax_idx);
            
            %current element
            Pnow=data.Points3D{it}; %Position of each vertex
            FCnow=FC{it,is,idef}; %Value of each face
            FANow = data.Deform.Area{it}; %Area of each face
            if optStruct.contactAreaLogic 
                Fcentroidnow=data.FaceCentroids{it}; %Position of each face centroids
            end
            
            %remove of unwanted elements
            FCnow(FCnow<optStruct.dataLimits(1))=NaN;
            FCnow(FCnow>optStruct.dataLimits(2))=NaN;
%             FCnow(abs(FCnow)>threshval(it,is,idef))=NaN; %Maximum value removal
            if optStruct.gapLogic == 1
                FCnow(data.FacePairInds==3)=NaN; %remove new mesh : too noisy 
            end
            AreaThreshold = prctile(FANow,99,'all');
            FCnow(FANow>AreaThreshold)=NaN; 
            
            %spatial filtering 
            if optStruct.smoothSpaceLogic
                FCnow=patchSmoothFaceMeasure(Fnow,Pnow,FCnow,optStruct);
            end
            
            if ~isempty(directionStringCell{is})
                Vnow=Vc{it,is};
                Dnow=Ds{it,is};
            end

            %contact area
            if optStruct.contactAreaLogic
                CAnow = find_contact_area(Fcentroidnow,Fcentroidlast,'lastmap_threshold',threshold_lastmap);
            end
            
            %------------------------------------------------------------------
            %%% ANIMATION STRUCTURE
            %------------------------------------------------------------------
            % FACE MEASURES (gpatch 1)
            animStruct.Handles{it}=[animStruct.Handles{it} hp(ax_idx) hp(ax_idx)];
            ii_anim=ii_anim+1;
            animStruct.Props{it}{ii_anim}='CData'; 
            animStruct.Set{it}{ii_anim}=FCnow; 
            ii_anim=ii_anim+1;
            animStruct.Props{it}{ii_anim}='Vertices';
            animStruct.Set{it}{ii_anim}=Pnow; 

            % FINGER INSIDE FILLING (gpatch 2)
            animStruct.Handles{it}=[animStruct.Handles{it} hp2(ax_idx)];
            ii_anim=ii_anim+1;
            animStruct.Props{it}{ii_anim}='Vertices';
            animStruct.Set{it}{ii_anim}=Pnow.*0.95+[0 -0.3 0.5]; 

            % CONTACT AREA DELIMITATION (plot3)
            if optStruct.contactAreaLogic
                animStruct.Handles{it}=[animStruct.Handles{it} hb(ax_idx) hb(ax_idx) hb(ax_idx)];
                ii_anim=ii_anim+1;
                animStruct.Props{it}{ii_anim}='XData'; 
                animStruct.Set{it}{ii_anim}=CAnow(:,1); 
                ii_anim=ii_anim+1;
                animStruct.Props{it}{ii_anim}='YData';
                animStruct.Set{it}{ii_anim}=CAnow(:,2);
                ii_anim=ii_anim+1;
                animStruct.Props{it}{ii_anim}='ZData'; 
                animStruct.Set{it}{ii_anim}=CAnow(:,3); 
            end
            
            % DIRECTIONS PRINCIPALES FACE MEASURE
            if ~isempty(directionStringCell{is})
                animStruct.Handles{it}=[animStruct.Handles{it} hq(ax_idx) hq(ax_idx) hq(ax_idx) hq(ax_idx) hq(ax_idx) hq(ax_idx)]; 
                ii_anim=ii_anim+1; 
                animStruct.Props{it}{ii_anim}='XData'; 
                animStruct.Set{it}{ii_anim}=Vnow(:,1); 
                ii_anim=ii_anim+1; 
                animStruct.Props{it}{ii_anim}='YData'; 
                animStruct.Set{it}{ii_anim}=Vnow(:,2); 
                ii_anim=ii_anim+1;
                animStruct.Props{it}{ii_anim}='ZData'; 
                animStruct.Set{it}{ii_anim}=Vnow(:,3); 
                ii_anim=ii_anim+1;
                animStruct.Props{it}{ii_anim}='UData'; 
                animStruct.Set{it}{ii_anim}=Dnow(:,1); 
                ii_anim=ii_anim+1; 
                animStruct.Props{it}{ii_anim}='VData'; 
                animStruct.Set{it}{ii_anim}=Dnow(:,2);
                ii_anim=ii_anim+1; 
                animStruct.Props{it}{ii_anim}='WData'; 
                animStruct.Set{it}{ii_anim}=Dnow(:,3);
            end
            ax_cur.XLim = xl; ax_cur.YLim = yl; ax_cur.ZLim = zl;
        end
    end
    
    % ROBOT RAW KINEMATICAL AND DYNAMICAL DATA  
    if optStruct.showRobotLogic 
        %current time in the transposed from the frame domain frequency
        tnow = data_robot.time(optStruct.frame_idx(it))-timestart;
        
        animStruct.Handles{it}=[animStruct.Handles{it},h_rule(1),h_rule(1),h_rule(2),h_rule(2)];
        
        ii_anim=ii_anim+1; 
        animStruct.Props{it}{ii_anim}='XData'; 
        animStruct.Set{it}{ii_anim}=[tnow tnow]; 
        ii_anim=ii_anim+1; 
        animStruct.Props{it}{ii_anim}='YData';
        animStruct.Set{it}{ii_anim}=[0 maxvalrule(1)]; 
        ii_anim=ii_anim+1; 
        animStruct.Props{it}{ii_anim}='XData'; 
        animStruct.Set{it}{ii_anim}=[tnow tnow]; 
        ii_anim=ii_anim+1; 
        animStruct.Props{it}{ii_anim}='YData';
        animStruct.Set{it}{ii_anim}=[0 maxvalrule(2)];
    end
    % RAW IMAGES FROM ONE CAMERA 
    if optStruct.showImgLogic 
        animStruct.Handles{it}=[animStruct.Handles{it},himg];
        
        ii_anim=ii_anim+1; 
        animStruct.Props{it}{ii_anim}='CData'; 
        animStruct.Set{it}{ii_anim}=data_img(:,:,it); 
    end
end

anim8(hf,animStruct);

%gui tools
addColorbarLimitsButton(hf);
addColormapButton(hf);
addEdgeColorButton(hf);
addFaceAlphaButton(hf);
addLightButton(hf);
addAmbientStrengthButton(hf);
addDiffuseStrengthButton(hf);
addSpecularStrengthButton(hf);
addQuiverFactorButton(hf);
addFaceLightingButton(hf);
set(hf,'visible',optStruct.showvisu);


if nargout == 1 
    varargout{1} = h_ax; 
end

end

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% % 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>