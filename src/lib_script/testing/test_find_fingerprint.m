

subject = 'S09'; 
material = 'coating'; 
trialnumber = '007'; 
camNbr = 2; 
trialname = subject+"_"+material+"_speckles_"+trialnumber;
file_cam = fullfile(baseDataPath,'rawdata',subject,"speckles",material,"vid",...
    sprintf('%s*_cam_%d*.mp4', trialname, camNbr));
S = dir(file_cam); 
vidraw = readvid(fullfile(S.folder, S.name));

%%
framerange = [1:150];
vid = vidraw(:,:,framerange); 
siz = size(vid); 
vidraw_1 = vid(:,:,1); 
% figure; imshow(vidraw_1(siz(1)/2:siz(1)/2+100,siz(2)/2:siz(2)/2+100),[]); 
y = prctile(vidraw_1(siz(1)/2:siz(1)/2+100,siz(2)/2:siz(2)/2+100),[1 99],'all');
y = double(y);
vid = satur(vid,'level',120); 

% mask = false(siz(1),siz(2)); 
% mask(100:siz(1)-100,100:siz(2)-100) = true; 
% 
% for ii = 1:siz(3)
%     vid_ii = vid(:,:,ii); 
%     vid_ii(~mask) = 0; 
%     vid(:,:,ii) = vid_ii; 
% end 
lowfilt = [25]; 
highfilt = [200,300]; 

vidfilt = cell(length(lowfilt),1); 
ivid=0; 
for ii = lowfilt
    ivid = ivid+1; 
    vidfilt{ivid} = filter_like_ben(vid,'gsbound',y,'paramfilt',[ii 300]);
    vidfilt{ivid}=satur(satur(vidfilt{ivid},'method','high','level',prctile(vidfilt{ivid}(:,:,end),99,'all')),...
                          'method','low','level',prctile(vidfilt{ivid}(:,:,end),1,'all')); 
end
% player(vidfilt{2}-vidfilt{1}); 
% %%
player(cat(1,[vidfilt{:}]))
%%
player(vidfilt{1});
%%
player(satur(satur(cat(1,vidfilt{:}),'method','high','level',prctile(vidfilt{end}(:,:,end),99,'all')),...
                          'method','low','level',prctile(vidfilt{end}(:,:,end),1,'all'))); 
%%
figure; imhist(vidraw_1)
%%
player(satur(vid,'level',120));

%%
% figure; imshow(vidraw(:,:,end)-vidraw(:,:,1)<20);
%%
player(satur(satur(vidfilt{end},'method','high','level',prctile(vidfilt{end}(:,:,end),95,'all')),...
                          'method','low','level',prctile(vidfilt{end}(:,:,end),5,'all'))); 



