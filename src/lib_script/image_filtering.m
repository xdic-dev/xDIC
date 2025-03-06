function imbdp = image_filtering(imarray,level)
imarray = satur(imarray,'level',level);

hsize = 11; sigma = 2; 
hpf=padarray(1,[(hsize-1)/2 (hsize-1)/2])-fspecial('gaussian',hsize,sigma);
lpf=fspecial('gaussian',hsize,sigma);

%replicate is used so that the image filtered doesn't suffer to much from
%border effect of the filtering process.
%
%im2single converts the intensity of the image to single precision 
imbdp=imfilter(imfilter(im2single(imarray),hpf,'replicate'),lpf,'replicate');

% figure; imshow(imbdp(:,:,1),[]); figure;imhist(imbdp(:,:,1)); 

% Remove extreme values from images
% imbdp = im_highlight(imbdp,refmask,'factor',10); 
y = prctile(imbdp(:,:,1),[45 55],'all'); 
% figure; imhist((imbdp(:,:,1)-y(1))/(y(2)-y(1)));

%rescale of the image to its original 8-bit precision  
imbdp = im2uint8((imbdp-y(1))/(y(2)-y(1))); 
% figure; imshow(imbdp(:,:,1),[]); figure;imhist(imbdp(:,:,1)); 

end