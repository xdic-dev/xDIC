function [C_smooth]=patchSmoothFaceMeasure_DD(varargin)

% function [C_smooth]=patchSmoothFaceMeasure(F,V,C,smoothPar)

%% Parse input
switch nargin 
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        smoothPar=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        smoothPar=varargin{4};
end

%% DD changed
smoothPar = smoothPar.smoothPar;
%% Get connectivity array

[connectivityStruct]=patchConnectivity(F,V);
faceFaceConnectivity=connectivityStruct.face.face;

%%

nDims=size(C,2); %Number of dimensions
logicValid=faceFaceConnectivity>0;
C_smooth=C;
C_smooth_step=C; 
for qIter=1:smoothPar.n 
    %Loop for all dimensions
    for qDim=1:1:nDims
        Xp=NaN(size(C,1),size(faceFaceConnectivity,2));
        Xp(logicValid)=C_smooth(faceFaceConnectivity(logicValid),qDim);
        Xp=mean(Xp,2,'omitnan');        
        C_smooth_step(:,qDim)=Xp;
    end
    C_smooth=((1-smoothPar.lambda).*C_smooth)+(smoothPar.lambda.*C_smooth_step);
end


