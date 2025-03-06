Nframe = length(DIC3DPPresults.Deform.J);
data = zeros(length(DIC3DPPresults.Deform.J{1}),length(DIC3DPPresults.Deform.J));
for ii = 1:Nframe
    data(:,ii) = DIC3DPPresults.Deform.J{ii}; 
end
x = data(4131,:); 

%filter param
[B,A] = butter(4,0.12); % low-pass filter
%
xcum = x; 
xcumf = filtfilt(B,A,xcum); 
%
xrate = [0,diff(x)]; 
xratef = filtfilt(B,A,xrate);
% double filter  
xratef2 = filtfilt(B,A,[0,diff(xcumf)]); 

figure; 
subplot(2,1,1); hold on; 
plot(xcum); 
plot(xcumf); 
plot(cumsum(xratef)+xcum(1)); 
legend('original','filtered direclty','reconstructed from filtered signal'); 

subplot(2,1,2); hold on; 
plot(xrate); 
plot(xratef); 
plot([0,diff(xcumf)]); 
plot(xratef2); 
legend('original','filtered direclty','reconstructed from filtered signal','double filtered'); 
