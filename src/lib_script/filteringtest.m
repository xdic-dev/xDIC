function im_modif = filteringtest(im)
siz = size(im); 
sigma = 1; 
level = 70; 
mask1 = imbinarize(im,level/255); %find the region of reflections
im_satur = satur(im,'level',level);
imSmooth = imgaussfilt3(im_satur, sigma,'FilterDomain','spatial','padding',5);
mask2 = (im-imSmooth)>20; 
mask = mask1|mask2;
% for level = [60,75]
% BW = false(siz);
% roi = false(siz);
% roi(:,:,1) = reference_save.roi.mask; 
im_modif = im_satur; 
% temproi = current_save(1).roi.mask; 
% p = zeros(siz(1)*siz(2),2,siz(3)); 
for ii = 1:siz(3)
%     if ii ~= 1
%     temproi = current_save(ii-1).roi.mask; 
%     end
    %BW(:,:,ii) = mask(:,:,ii) ;%| boundarymask(mask(:,:,ii),4);
    tempmask = mask(:,:,ii);%|boundarymask(mask(:,:,ii),8);
    tempim = im_satur(:,:,ii);
    tempsmooth = imSmooth(:,:,ii); 
%     tempsmooth(~temproi)=0; 
%     tempim(~temproi)=0;
    tempim(tempmask)=tempsmooth(tempmask); %replace the reflections by smoothed pixels 
    im_modif(:,:,ii)=tempim;
%     [y,x] = find(tempmask);
%     p(1:length(y),:,ii) = [x,y]; 
end
end