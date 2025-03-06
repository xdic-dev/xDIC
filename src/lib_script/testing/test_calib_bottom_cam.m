%%% TEST CALIBRATION JAI 
subject = 'S09'; 
trialnumber = '013'; 
material = 'coating'; 
camNbr = 0; 
trialname = subject+"_"+material+"_speckles_"+trialnumber;

file_cam = fullfile(baseDataPath,'rawdata',subject,"speckles",material,"vid",...
    sprintf('%s*_cam_%d*.mp4',trialname,camNbr));

frameStart = 1;
frameEnd = 150; 
framejump = 1; 
S = dir(file_cam);
if isempty(S)
    error('Error: Video file not found.');
end
cam_raw_1 = readvid(fullfile(S(1).folder, S(1).name),frameEnd-frameStart+1,frameStart,framejump);
cam_raw_2 = readvid(fullfile(S(2).folder, S(2).name),frameEnd-frameStart+1,frameStart,framejump);


file_cam = fullfile(baseDataPath,'rawdata',subject,"calib",'2',...
    sprintf('*_%d.tiff',camNbr));
S = dir(file_cam);
cam_calib_1 = imread(fullfile(S(1).folder,S(1).name));
file_cam = fullfile(baseDataPath,'rawdata',subject,"calib",'1',...
    sprintf('*_%d.tiff',camNbr));
S = dir(file_cam);
cam_calib_2 = imread(fullfile(S(end).folder,S(end).name));
%% figure;
figure; 
imshow(cat(2,cam_raw_1(:,:,end),cam_raw_2(:,:,end)));
%%
player(permute(cam_raw_1,[2 1 3])); 
%%
smallImage = imresize(cam_calib_1, 1/2, 'bilinear');
smallImage = flip(smallImage',2);
figure; 
imshow(smallImage); hold on; 
H = drawcircle; hold on; 
pause; 
c = H.Center; 
r = H.Radius;
drawpoint(gca,'Position',[c(1),c(2)]); hold on; 
%
H = drawline; hold on; 
pause;
%
p1 = H.Position(1,:)'; 
p2 = -H.Position(2,:)'; 
alpha = atan2((p2(2)-p1(2)),(p2(1)-p1(1)))*180/pi
%%
figure; 
imshow(cam_calib_2); hold on; 
H = drawline; hold on; 
pause;
p1 = H.Position(:,1)'; 
p2 = H.Position(:,2)'; 
d = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
nbrMM = 23; %%mm from left to right
res = d/nbrMM



