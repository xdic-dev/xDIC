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

% %%
% 
% nDims=size(C,2); %Number of dimensions
% logicValid=faceFaceConnectivity>0;
% C_smooth=C;
% C_smooth_step=C; 
% for qIter=1:smoothPar.n 
%     %Loop for all dimensions
%     for qDim=1:1:nDims
%         Xp=NaN(size(C,1),size(faceFaceConnectivity,2));
%         Xp(logicValid)=C_smooth(faceFaceConnectivity(logicValid),qDim);
%         Xp=mean(Xp,2,'omitnan');        
%         C_smooth_step(:,qDim)=Xp;
%     end
%     C_smooth=((1-smoothPar.lambda).*C_smooth)+(smoothPar.lambda.*C_smooth_step);
% end

N = smoothPar.n; 
sigma = smoothPar.sigma; 

Co = faceFaceConnectivity; 
sizCo = size(Co); 
sizC = size(C);
W = gaussmf(0:N,[sigma 0]); %W = W/sum(W); %W = flip(W/sum(W)); 

% figure; stem(W); 

%stage 0 : only the first face 
idx_kk=1:sizC; 
Cs_adj = C(idx_kk)*W(1); %stage 0 %zeros(sizC); 

%stage 1 : 1st-adjacent faces are used.
if N >= 1
    Xpset = NaN(sizCo(1),1);
    Xp=NaN(sizCo(1),sizCo(2));
    
    % Get centered indices of faces to evaluate
    idx_kk=1:sizC; 
    
    % Get centered face value 
    Xpset(idx_kk) = C(idx_kk); 
    
    %current face connectivity matrix
    Co_kk=Co(idx_kk,:); %initial one at this stage
    
    %current valid face connectivity matrix
    logicValid_kk=Co_kk>0; %non NaN, takes only existing faces
    
    % Retreive face measure of adjacent faces
    Xp(logicValid_kk)=C(Co_kk(logicValid_kk),1);
    
    % Mean of current and adjacent faces measures 
    Xp=mean([Xpset,Xp],2,'omitnan');
    
    Cs_adj = Cs_adj+Xp*W(2);
end

%stage 2 : N-th adjacent faces are used. 
% initialized parameters : Co_kk, logicValid_kk
for kk = 2:N
    Xpm = zeros(sizCo(1),sizCo(2)); 
    for ii = 1:sizCo(2)
        Xpset = NaN(sizCo(1),1);
        Xp=NaN(sizCo(1),sizCo(2));
        
        % Get centered indices of faces to evaluate
        idx_kk = Co_kk(logicValid_kk(:,ii),ii); 
        
        % Get centered face value 
        Xpset(idx_kk) = C(idx_kk); 
        
        %current face connectivity matrix
        Co_kk = Co(idx_kk,:); 
        
        %current valid face connectivity matrix
        logicValid_kk=Co_kk>0;
        
        % Retreive face measure from adjacent faces
        Xp(logicValid_kk)=C(Co_kk(logicValid_kk),1); 
        
        size(C(idx_kk))
        size(Xp)
        % Mean of current and adjacent faces measures 
        Xp=mean([Xpset,Xp],2,'omitnan');
        
        Xpm(:,ii)=Xp; 
    end
    Xpm=mean(Xpm,2,'omitnan');
    Cs_adj=Cs_adj+Xpm*W(kk+1); %stage 2
end

% Weighted sum over all contributions of adjacent faces. 
C_smooth = Cs_adj/sum(W(1:N));
end






















