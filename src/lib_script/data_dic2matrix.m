function [FCmat_cum,FCmat_rate] = data_dic2matrix(data, faceMeasureCell_is)
nFrames = length(data.Disp.DispVec); 
if strcmp(faceMeasureCell_is,'DispX')
    F=data.Faces;
    dispNow=data.Disp.DispVec{1}(:,1); % point measure
    z = zeros(nFrames,length(mean(dispNow(F),2)));
    FCmat = z;
    for it = 1:nFrames
        dispNow=data.Disp.DispVec{it}(:,1); % point measure
        FCmat(it,:)=mean(dispNow(F),2);
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
elseif strcmp(faceMeasureCell_is,'DispY')
    F=data.Faces;
    dispNow=data.Disp.DispVec{1}(:,2); % point measure
    z = zeros(nFrames,length(mean(dispNow(F),2)));
    FCmat = z;
    for it = 1:nFrames
        dispNow=data.Disp.DispVec{it}(:,2); % point measure
        FCmat(it,:)=mean(dispNow(F),2);
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
elseif strcmp(faceMeasureCell_is,'DispXY')
    F=data.Faces;
    dispNow=data.Disp.DispVec{1}(:,1); % point measure
    z = zeros(nFrames,length(mean(dispNow(F),2)));
    FCmat = z;
    for it = 1:nFrames
        dispNow=data.Disp.DispVec{it}(:,1); % point measure
        dispNow=dispNow+data.Disp.DispVec{it}(:,2); % point measure
        FCmat(it,:)=mean(dispNow(F),2);
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
elseif strcmp(faceMeasureCell_is,'DispZ')
    F=data.Faces;
    dispNow=data.Disp.DispVec{1}(:,3); % point measure
    z = zeros(nFrames,length(mean(dispNow(F),2)));
    FCmat = z;
    for it = 1:nFrames
        dispNow=data.Disp.DispVec{it}(:,3); % point measure
        FCmat(it,:)=mean(dispNow(F),2);
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';%-FCmat(24,:);
elseif strcmp(faceMeasureCell_is,'DispMgn')
    F=data.Faces;
    dispNow=data.Disp.DispMgn{1}; % point measure
    z = zeros(nFrames,length(mean(dispNow(F),2)));
    FCmat = z;
    for it = 1:nFrames
        dispNow=data.Disp.DispMgn{it}; % point measure
        FCmat(it,:)=mean(dispNow(F),2);
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
elseif strcmp(faceMeasureCell_is,'corrComb')
    z = zeros(nFrames,length(data.(faceMeasureCell_is){1}));
    FCmat = z;
    for it = 1:nFrames
        FCmat(it,:)=data.(faceMeasureCell_is){it};
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
elseif strcmp(faceMeasureCell_is,'FaceColors')
    z = zeros(nFrames,length(data.(faceMeasureCell_is)));
    FCmat = z;
    speckles_pattern = 255*(randi(2,length(data.(faceMeasureCell_is)),1)-1);
    for it = 1:nFrames
        FCmat(it,:)=speckles_pattern; %data.FaceColors; %speckles_pattern;
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
elseif isfield(data,faceMeasureCell_is) 
    z = zeros(nFrames,length(data.(faceMeasureCell_is)));
    FCmat = z;
    for it = 1:nFrames
        FCmat(it,:)=data.(faceMeasureCell_is);
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
else
    z = zeros(nFrames,length(data.Deform.(faceMeasureCell_is){1}));
    FCmat = z;
    for it = 1:nFrames
        FCmat(it,:)=data.Deform.(faceMeasureCell_is){it};
    end
    FCmat_rate = [z(1,:);diff(FCmat,1,1)]';
    FCmat_cum = FCmat';
end
end