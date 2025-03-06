function varargout = filter3Ddeform_time(varargin)
p = inputParser;
p.addRequired('data');
p.addParameter('freqFilt',10);
% p.addParameter('freqFilt_rate',10);
p.addParameter('freqAcq',50);
p.addParameter('F',[]);
p.addParameter('P',[]);
p.addParameter('smoothPar',[]);
p.addParameter('smoothSpace',0);
p.addParameter('deformCol',1);


p.parse(varargin{:});
data = p.Results.data;
freqFilt = p.Results.freqFilt;
% freqFilt_rate = p.Results.freqFilt_rate;
freqAcq = p.Results.freqAcq;
F = p.Results.F;
P = p.Results.P;
smoothPar = p.Results.smoothPar;
smoothSpace = p.Results.smoothSpace;
deformCol = p.Results.deformCol;

Nframe = length(data); 
Npoints = length(data{1}); 
Ncols = size(data{1},2); 
d = zeros(Nframe,Npoints); 
dcum_cell = cell(Nframe,1); 
% drate_cell = cell(Nframe,1); 
n = NaN(Npoints,1); 
for ii = 1:Nframe
    d(ii,:) = data{ii}(:,deformCol); 
end

%% Fill NaN - to be analyzed 
d = fillmissing(d,'linear',1,'EndValues','nearest'); 

%% Retreive non NaN values
mask = ~any(isnan(d),1); %all data point with no NaN at any point in time

%% Filtering 
npad = 3; % decrease border effect of the filter
% CUM FILTER
x = [repmat(d(1,mask),[npad,1]);...
    d(:,mask);...
    repmat(d(end,mask),[npad,1])]; 
[B,A] = butter(4,freqFilt/(freqAcq/2)); % low-pass filter 
data_cum_f = filtfilt(B,A,x);

% % RATE FILTER
% x = diff([data_cum_f(1,:);data_cum_f],1);
% [B,A] = butter(4,freqFilt_rate/(freqAcq/2)); % low-pass filter
% data_rate_f = filtfilt(B,A,x); 

%% Output for plot : 
% mat 2 cell conversion 
smoothPar.smoothPar = smoothPar; 
v = 1:Ncols; v(deformCol)=[]; 
for ii = 1:Nframe
    dcum_cell{ii} = n; %fill with NaN for the correct vector size 
%     drate_cell{ii} = n; %fill with NaN for the correct vector size 
    if smoothSpace == 1
        data_cum_f(ii+npad,:) = patchSmoothFaceMeasure(F,P{ii},data_cum_f(ii+npad,:)',smoothPar);
%         data_rate_f(ii+npad,:) = patchSmoothFaceMeasure(F,P{ii},data_rate_f(ii+npad,:)',smoothPar);
    end
    dcum_cell{ii}(mask,deformCol) = data_cum_f(ii+npad,:); 
    dcum_cell{ii}(mask,v) = data{ii}(mask,v);
%     drate_cell{ii}(mask,deformCol) = data_rate_f(ii+npad,:);
%     drate_cell{ii}(mask,v) = data{ii}(mask,v);
end
% first frame at 0 
dcum_cell{1}(mask) = zeros(size(d(1,mask)));
% drate_cell{1}(mask) = zeros(size(d(1,mask)));

%% output 
varargout{1} = dcum_cell; 
% varargout{2} = drate_cell; 
end
