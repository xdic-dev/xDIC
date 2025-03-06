function varargout = anim8_DIC3DPP_faceMeasure(DIC3DPPresults,faceMeasureString,RBMlogic,varargin)
%% function for plotting 3D-DIC results of face measures in STEP4.
% plotting 3D surfaces from camera pairs, animation changing
% with time, and the faces colored according to faceMeasureString
% this function is called in plotMultiDICPairResults
%
% Options:
% anim8_DIC3DPP_faceMeasure(D(DIC3DPPresults,faceMeasureString,RBMlogic)
% anim8_DIC3DPP_faceMeasure((DIC3DPPresults,faceMeasureString,RBMlogic,optStruct)
%
% Inputs:
% * DIC3DAllPairsResults
% * faceMeasureString: can be any of the following:
%   'J','Emgn','emgn','Epc1','Epc2','epc1','epc2','dispMgn','dispX','dispY','dispZ','FaceColors','FaceIsoInd','pairInd','Lamda1','Lamda2'
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
if ~isfield(optStruct,'smoothLogic')
    optStruct.smoothLogic=0;
end
if ~isfield(optStruct,'dataLimits')
    optStruct.dataLimits=[-inf inf];
end
if ~isfield(optStruct,'zDirection') % 1 or -1
    optStruct.zDirection=1;
end

if ~isfield(optStruct,'maxCorrCoeff')
    optStruct.maxCorrCoeff=[];
end
if ~isfield(optStruct,'gap')
    optStruct.gap=0;
end

%%
nFrames=numel(DIC3DPPresults.Points3D);
[xl,yl,zl]=axesLimits(DIC3DPPresults.Points3D);

%% Assign the right face measure into FC
if iscell(faceMeasureString)
    faceMeasureCell=faceMeasureString;
    nStrains=numel(faceMeasureString);
    FC=cell(nFrames,nStrains);
    
elseif ischar(faceMeasureString)
    nStrains=1;
    faceMeasureCell=cell(1);
    faceMeasureCell{1}=faceMeasureString;
    FC=cell(nFrames,1);
else
    error('wrong face measure (second input variable)');
end

for is=1:nStrains
switch faceMeasureCell{is}
    case 'FaceColors'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.(faceMeasureString);
            FC{it,is}(~(DIC3DPPresults.FaceCentroids{it}(:,3) > 0.5))=NaN;
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'FacePairInds'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.FacePairInds;
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'J'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it}; % face color (strain)
            if ~isempty(optStruct.maxabs_tot)
                if it ~= 1
                    FClast = FC{it-1,is}; 
                    mask1 = abs(FC{it,is}-FClast)>optStruct.maxabs_rate; 
                    mask2 = abs(FC{it,is})>optStruct.maxabs_tot;
                    FC{it,is}(mask1|mask2,:)=NaN;
                end
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'Lamda1'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxabs_tot)
                if it ~= 1
                    FClast = FC{it-1,is}; 
                    mask1 = abs(FC{it,is}-FClast)>optStruct.maxabs_rate; 
                    mask2 = abs(FC{it,is}-1)>optStruct.maxabs_tot;
                    FC{it,is}(mask1|mask2,:)=NaN;
                end
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'Lamda2'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxabs_tot)
                if it ~= 1
                    FClast = FC{it-1,is}; 
                    mask1 = abs(FC{it,is}-FClast)>optStruct.maxabs_rate; 
                    mask2 = abs(FC{it,is}-1)>optStruct.maxabs_tot;
                    FC{it,is}(mask1|mask2,:)=NaN;
                end
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'Emgn'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'emgn'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'Eeq'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'eeq'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'EShearMax'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it}; % face color (strain)
            if ~isempty(optStruct.maxabs_tot)
                if it ~= 1
                    FClast = FC{it-1,is}; 
                    mask1 = abs(FC{it,is}-FClast)>optStruct.maxabs_rate; 
                    mask2 = abs(FC{it,is})>optStruct.maxabs_tot;
                    FC{it,is}(mask1|mask2,:)=NaN;
                end
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'eShearMax'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'Epc1'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxabs_tot)
                if it ~= 1
                    FClast = FC{it-1,is}; 
                    mask1 = abs(FC{it,is}-FClast)>optStruct.maxabs_rate; 
                    mask2 = abs(FC{it,is})>optStruct.maxabs_tot;
                    FC{it,is}(mask1|mask2,:)=NaN;
                end
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'Epc2'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxabs_tot)
                if it ~= 1
                    FClast = FC{it-1,is}; 
                    mask1 = abs(FC{it,is}-FClast)>optStruct.maxabs_rate; 
                    mask2 = abs(FC{it,is})>optStruct.maxabs_tot;
                    FC{it,is}(mask1|mask2,:)=NaN;
                end
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'epc1'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'epc2'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'DispMgn'
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            dispNow=DIC3DPPresults.Disp.DispMgn{it}; % point measure
            FC{it,is}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'DispX'
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            dispNow=DIC3DPPresults.Disp.DispVec{it}(:,1); % point measure
            FC{it,is}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'DispY'
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            dispNow=DIC3DPPresults.Disp.DispVec{it}(:,2); % point measure
            FC{it,is}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'DispZ'
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
%             dispNow=DIC3DPPresults.Disp.DispVec{it}(:,3)-DIC3DPPresults.Disp.DispVec{17}(:,3); % point measure
            dispNow=DIC3DPPresults.Points3D{it}(:,3); % point measure
            FC{it,is}=mean(dispNow(Fnow),2); % turn into face measure
            if isfield(optStruct,'zlim') % 'none' or 'k'
                FC{it,is}(FC{it,is}<optStruct.zlim,:) = NaN;
            end
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'FaceIsoInd'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
    case 'FaceCorrComb'
        for it=1:nFrames
            FC{it,is}=DIC3DPPresults.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        
    otherwise
        error('unexpected face measure string. plots not created');
end     
end

%% Plot
%figure
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

nrows = 1; 
nPropsSupp = 0; 
data_robot = optStruct.data_robot;
showrobot = ~isempty(data_robot);
if showrobot 
    nrows = 2; 
    nPropsSupp = nrows*2; %two data per 2Dplot : X and Y 
end

mergevec = mat2cell((repelem([1;nStrains+2],1,nStrains)+(0:nStrains-1))',repelem(1,nStrains),2)';
h_ax = subplot_ax(nrows,nStrains+(nrows-1),'merge',mergevec,'shownum',0); hold(h_ax,'on'); 
set(h_ax,'FontSize',15);

for is=1:nStrains
    %current subplot
    ax_idx = is; 
    set(hf,'CurrentAxes',h_ax(ax_idx)); 
    %axis
    axisGeom(h_ax(ax_idx));
    %view
    h_ax(ax_idx).CameraUpVector=[0 0 optStruct.zDirection];
    %colorbar
%     caxis(h_ax(ax_idx),optStruct.colorBarLimits);
    colorbar(h_ax(ax_idx));
    %title
    title(h_ax(ax_idx),optStruct.supTitleString);

    it=1;
    Fnow=DIC3DPPresults.Faces;
    Fcentroidlast=DIC3DPPresults.FaceCentroids{end};
    Fcentroidnow=DIC3DPPresults.FaceCentroids{it};
    Pnow=DIC3DPPresults.Points3D{it};
    FCnow=FC{it};
    if optStruct.smoothLogic
        [FCnow]=patchSmoothFaceMeasure(Fnow,Pnow,FCnow,optStruct);
    end
    FCnow(FCnow<optStruct.dataLimits(1))=NaN;
    FCnow(FCnow>optStruct.dataLimits(2))=NaN;    
    if optStruct.gap == 1
        FCnow(DIC3DPPresults.FacePairInds==3)=NaN; %remove new mesh : too noisy 
    end

    threshold = -1; 
    %refine threshold from last map 
    idx_IN = Fcentroidlast(:,3) < threshold; 
    DM = [Fcentroidlast(idx_IN,1),Fcentroidlast(idx_IN,2),ones(length(Fcentroidlast(idx_IN,1)),1)]; % Design Matrix
    param = DM\Fcentroidlast(idx_IN,3);                   % Estimate Parameters
    B = param; 
    planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.13);
    plane_id = planefit(Fcentroidlast(:,1),Fcentroidlast(:,2)) > Fcentroidlast(:,3); 
    threshold = mean(Fcentroidlast(plane_id,3));

    %fit plane with point below threshold
    IN_threshold = 300; 
    idx_IN = Fcentroidnow(:,3) < threshold+0.5; 
    DM = [Fcentroidnow(idx_IN,1),Fcentroidnow(idx_IN,2),ones(length(Fcentroidnow(idx_IN,1)),1)]; % Design Matrix
    param = DM\Fcentroidnow(idx_IN,3);                   % Estimate Parameters
    B = param; 
    planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.20);

    %find the boundary 
    plane_id = planefit(Fcentroidnow(:,1),Fcentroidnow(:,2)) > Fcentroidnow(:,3); 
    if threshold+0.5>B(3)
        CAnow = findboundary(Fcentroidnow(plane_id,:));
    else
        CAnow = [NaN,NaN,NaN];
    end

    %plot
    hp(ax_idx)=gpatch(Fnow,Pnow,FCnow,optStruct.lineColor,optStruct.FaceAlpha); hold on
	hp2(ax_idx)=gpatch(Fnow,Pnow.*[0.95 1 1]+[0 -0.1 0.05]); hold on
    hb(ax_idx)=plot3(CAnow(:,1),CAnow(:,2),CAnow(:,3),'k-','LineWidth',5); hold on
    %
    h_ax(ax_idx).XLim = xl; h_ax(ax_idx).YLim = yl; h_ax(ax_idx).ZLim = zl;
end
%robot plot
if showrobot
    %displacement 
    ax_idx = 1*(nStrains+1); 
    set(hf,'CurrentAxes',h_ax(ax_idx)); 

    timestart = data_robot.time(optStruct.frame_idx(1));
    timeend = data_robot.time(optStruct.frame_idx(end));
    timeframe = optStruct.frame_idx(1):optStruct.frame_idx(end);
    dispX = data_robot.position_robot(timeframe,1)-data_robot.position_robot(optStruct.frame_idx(1),1);
    dispY = data_robot.position_robot(timeframe,2)-data_robot.position_robot(optStruct.frame_idx(1),2);
    disp = sqrt(dispX.^2+dispY.^2);
    maxvalrule(1) = ceil(max(disp)); 
    %
    plot(h_ax(ax_idx),data_robot.time(timeframe)-timestart,disp,'k-','LineWidth',3);
    %
    h_rule(1)=plot(h_ax(ax_idx),[0 0],[0 maxvalrule(1)],'r-',...
        'LineWidth',5);
    %
    h_ax(ax_idx).XLabel.String = 'Time (s)'; h_ax(ax_idx).YLabel.String = 'disp (mm)';
    h_ax(ax_idx).XLim = [0 timeend-timestart]; h_ax(ax_idx).YLim = [0 maxvalrule(1)];

    %force 
    ax_idx = 2*(nStrains+1); 
    set(hf,'CurrentAxes',h_ax(ax_idx)); 

    timestart = data_robot.time(optStruct.frame_idx(1));
    timeend = data_robot.time(optStruct.frame_idx(end));
    maxvalrule(2) = ceil(max(data_robot.tf(optStruct.frame_idx))); 
    if maxvalrule(2)<ceil(max(data_robot.nf(optStruct.frame_idx)))
        maxvalrule(2) = ceil(max(data_robot.nf(optStruct.frame_idx)));
    end
    timeframe = optStruct.frame_idx(1):optStruct.frame_idx(end);
    %
    plot(h_ax(ax_idx),data_robot.time(timeframe)-timestart,data_robot.nf(timeframe),'k-','LineWidth',3);
    plot(h_ax(ax_idx),data_robot.time(timeframe)-timestart,data_robot.tf(timeframe),'k-','LineWidth',3,'Color',[1 1 1].*0.5);
    %
    h_rule(2)=plot(h_ax(ax_idx),[0 0],[0 maxvalrule(2)],'r-',...
        'LineWidth',5);
    %
    h_ax(ax_idx).XLabel.String = 'Time (s)'; h_ax(ax_idx).YLabel.String = 'Force (N)';
    h_ax(ax_idx).XLim = [0 timeend-timestart]; h_ax(ax_idx).YLim = [0 maxvalrule(2)];
    % 
    legend(h_ax(ax_idx),{'NF','TF'}); 
end

%% fill in the animstuct for each frame
animStruct.Time=1:nFrames;
animStruct.Handles=cell(1,nFrames);
animStruct.Props=cell(1,nFrames);
animStruct.Set=cell(1,nFrames);

for it=1:nFrames
    nProps = 7; 
    animStruct.Handles{it}=[];
    animStruct.Props{it}=cell(1,nProps*nStrains+nPropsSupp);
    animStruct.Set{it}=cell(1,nProps*nStrains+nPropsSupp);
    
    for is=1:nStrains
        ax = h_ax(is);
        
        Pnow=DIC3DPPresults.Points3D{it};
        Fnow=DIC3DPPresults.Faces;
        Fcentroidnow=DIC3DPPresults.FaceCentroids{it};
        FCnow=FC{it,is};
        if optStruct.smoothLogic
            [FCnow]=patchSmoothFaceMeasure(Fnow,Pnow,FCnow,optStruct);
        end

        FCnow(FCnow<optStruct.dataLimits(1))=NaN;
        FCnow(FCnow>optStruct.dataLimits(2))=NaN; 
        if optStruct.gap == 1
            FCnow(DIC3DPPresults.FacePairInds==3)=NaN; %remove new mesh : too noisy 
        end

        %fit plane with point below threshold
        idx_IN = Fcentroidnow(:,3) < threshold+0.5; 
        DM = [Fcentroidnow(idx_IN,1),Fcentroidnow(idx_IN,2),ones(length(Fcentroidnow(idx_IN,1)),1)]; % Design Matrix
        param = DM\Fcentroidnow(idx_IN,3);                   % Estimate Parameters
        B = param;
        planefit = @(x,y) B(2)*y + (B(3)+0.00); %B(1)*x +  

        %find the boundary 
        plane_id = planefit(Fcentroidnow(:,1),Fcentroidnow(:,2)) > Fcentroidnow(:,3); 
        if threshold+0.5>B(3)
            CAnow = findboundary(Fcentroidnow(plane_id,:));
        else
            CAnow = [NaN,NaN,NaN];
        end

        animStruct.Handles{it}=[animStruct.Handles{it} hp(is) hp(is) hp2(is) hb(is) hb(is) hb(is)]; %Handles of objects to animate
%         animStruct.Handles{it}=[animStruct.Handles{it} hp(is) hp(is) hb(is) hb(is) hb(is)]; %Handles of objects to animate
        ii = 1; 
        animStruct.Props{it}{ii+nProps*(is-1)}='CData'; ii=ii+1; 
        animStruct.Props{it}{ii+nProps*(is-1)}='Vertices'; ii=ii+1;%Properties of objects to animate
        animStruct.Props{it}{ii+nProps*(is-1)}='Vertices'; ii=ii+1;%Properties of objects to animate
        animStruct.Props{it}{ii+nProps*(is-1)}='XData'; ii=ii+1;%Properties of objects to animate
        animStruct.Props{it}{ii+nProps*(is-1)}='YData'; ii=ii+1;%Properties of objects to animate
        animStruct.Props{it}{ii+nProps*(is-1)}='ZData'; ii=ii+1;%Properties of objects to animate

        ii = 1; 
        animStruct.Set{it}{ii+nProps*(is-1)}=FCnow;ii=ii+1;
        animStruct.Set{it}{ii+nProps*(is-1)}=Pnow; ii=ii+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+nProps*(is-1)}=Pnow.*[0.95 0.95 0.95]; ii=ii+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+nProps*(is-1)}=CAnow(:,1); ii=ii+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+nProps*(is-1)}=CAnow(:,2); ii=ii+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+nProps*(is-1)}=CAnow(:,3); ii=ii+1;%Property values for to set in order to animate
       
        ax.XLim = xl; ax.YLim = yl; ax.ZLim = zl;   
    end
    
    if showrobot
        jj = 1; ii = ii-1; 
        animStruct.Handles{it}=[animStruct.Handles{it},h_rule(1),h_rule(1),h_rule(2),h_rule(2)];
        animStruct.Props{it}{ii+jj+nProps*(nStrains-1)}='XData'; jj=jj+1;%Properties of objects to animate
        animStruct.Props{it}{ii+jj+nProps*(nStrains-1)}='YData'; jj=jj+1;%Properties of objects to animate
        animStruct.Props{it}{ii+jj+nProps*(nStrains-1)}='XData'; jj=jj+1;%Properties of objects to animate
        animStruct.Props{it}{ii+jj+nProps*(nStrains-1)}='YData'; jj=jj+1;%Properties of objects to animate
        
        jj=1;
        timenow = data_robot.time(optStruct.frame_idx(it))-timestart;
        animStruct.Set{it}{ii+jj+nProps*(nStrains-1)}=[timenow timenow]; jj=jj+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+jj+nProps*(nStrains-1)}=[0 maxvalrule(1)]; jj=jj+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+jj+nProps*(nStrains-1)}=[timenow timenow]; jj=jj+1;%Property values for to set in order to animate
        animStruct.Set{it}{ii+jj+nProps*(nStrains-1)}=[0 maxvalrule(2)]; jj=jj+1;%Property values for to set in order to animate

        ax_idx = 1*(nStrains+1); 
        h_ax(ax_idx).XLim = [0 timeend-timestart]; h_ax(ax_idx).YLim = [0 maxvalrule(1)];
        ax_idx = 2*(nStrains+1); 
        h_ax(ax_idx).XLim = [0 timeend-timestart]; h_ax(ax_idx).YLim = [0 maxvalrule(2)];
    end
end

% %%% DD 
% % Change edge color for all existing images accordingly
% hpatches = findobj(hf,'type','patch');
% numPatches = size(hpatches,1);
% for i = 1:numPatches
%     hpatches(i).EdgeColor = 'none';
% end
% %%%
anim8(hf,animStruct);

addColorbarLimitsButton(hf);
addColormapButton(hf);
addEdgeColorButton(hf);
addFaceAlphaButton(hf);
addLightButton(hf);
addAmbientStrengthButton(hf);
addDiffuseStrengthButton(hf);
addSpecularStrengthButton(hf);
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
%
% Modified by Rana Odabas 2018
%
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>