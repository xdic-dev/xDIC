
cell = vidraw{2}; 
for ii = 1:2
    subject = 'S09'; 
    material = 'coating'; 
    trialnumber = '007'; 
    camNbr = ii; 
    trialname = subject+"_"+material+"_speckles_"+trialnumber;
    file_cam = fullfile(baseDataPath,'rawdata',subject,"speckles",material,"vid",...
        sprintf('%s*_cam_%d*.mp4', trialname, camNbr));
    S = dir(file_cam); 
    vidraw{1} = readvid(fullfile(S.folder, S.name));
end
vidraw1 = vidraw{1};
vidraw2 = vidraw{2};

%% Filtering
%mask around the center of the image to take out the greyscale value of a
%patch containing speckles. 
siz = size(vidraw1);
mask = false(siz(1),siz(2));
mask(siz(1)/2-50:siz(1)/2+50,siz(2)/2-50:siz(2)/2+50) = true;

vidfilt1 = filter_like_ben(vidraw1(:,:,1:150),'mask',mask,'paramfilt',[20,300]);
vidfilt2 = filter_like_ben(vidraw2(:,:,1:150),'mask',mask,'paramfilt',[20,300]);

% MASK BORDER
% siz = size(im1); 
% maskborder = false(size(im1));
% maskborder(100:size(im1,1)-100,100:size(im1,2)-100,:) = true; 
% im1(~maskborder) = 0; im1 = reshape(im1,siz); 
% im2(~maskborder) = 0; im2 = reshape(im2,siz); 

%%
% FIND POSITION OF SUBSET FROM MAX CORRELATION VALUE  
useFiltImgLogic = 0; 
for ii = 1:2
if ii == 1
    im2 = vidfilt2(:,:,1); % ref seed subset
    im1 = vidfilt1(:,:,1); % def image
else
    im2 = vidraw2(:,:,1); % ref seed subset
    im1 = vidraw1(:,:,1); % def image
end



useLastRectangleLogic = ii == 2; 
if ~useLastRectangleLogic
    % DRAW RECTANGLE 
    fig = figure; 
    title('Draw subset to be found') 
    imshow(im2); 
    hrec = drawrectangle(gca);
    initial_ROI = get(hrec,'Position'); initial_ROI = round(initial_ROI); 
    xstart=initial_ROI(1);ystart=initial_ROI(2);rec_w=initial_ROI(3);rec_h=initial_ROI(4);
    close(fig)
end

% FIND BEST CORRELATION BETWEEN SUBSET AND INITIAL IMAGE. 
im2_subset = im2(ystart:ystart+rec_h,xstart:xstart+rec_w);
siz = size(im2_subset); 

c = normxcorr2(im2_subset,im1); 
stdCorrCoeff = std(c(:)); 
[maxCorrCoeff,id] = max(c(:)); 
[ypeak,xpeak] = find(c==maxCorrCoeff);
cz = imresize(c(round(ypeak-rec_h/2:ypeak+rec_h/2),round(xpeak-rec_w/2:xpeak+rec_w/2)),[NaN,size(c,2)]);
figure; 
subplot(2,2,1); hold on; imshow([c;cz],[]); colormap('jet');  axis equal; shading faceted
subplot(2,2,2); hold on; plot(1:length(c(:)),c(:),'.'); plot(id,maxCorrCoeff,'ro'); ylabel('CorrCoeff');
if maxCorrCoeff > 0.3
    yoffSet = ypeak-siz(1);
    xoffSet = xpeak-siz(2);
    subplot(2,2,3); imshow([im2;im1]); hold on; 
    title(sprintf('[%d,%d]->[%d,%d]',ystart,xstart,ypeak,xpeak));  hold on; 
    drawrectangle(gca,'Position',initial_ROI, ...
        'FaceAlpha',0);
    drawrectangle(gca,'Position',[xoffSet,yoffSet+size(im2,1),siz(2),siz(1)], ...
        'FaceAlpha',0);
    im1_subset_found = im1(yoffSet:yoffSet+siz(1)-1,xoffSet:xoffSet+siz(2)-1);
    subplot(2,2,4); 
    imshow([im2_subset;im1_subset_found]); hold on; 
    title(sprintf('MaxCorrCoeff : %1.3f\n STD : %1.3f / Ratio : %1.3f',maxCorrCoeff,stdCorrCoeff,maxCorrCoeff/stdCorrCoeff));  
else 
    fprintf('Seed not found try again\n'); 
end
end