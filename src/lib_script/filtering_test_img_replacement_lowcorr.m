%filename
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 
%target trial info 
subject = 'S01'; 
fingerpattern = 'speckles'; 
material = 'coating'; %'coating'; %'glass';
trialnumber = '001'; 
phase = 'loading'; %'loading'; %'slide1';

filename = 'ncorr1.mat';
file_ncorr = fullfile(baseResultPath,'analysis',subject,material,trialnumber,phase,'test_30radius',filename);
%load results
r = load(file_ncorr); 
%%
Nframe = length(r.current_save);
siz = [size(r.current_save(1).gs),Nframe]; 
z = zeros(siz); 
im=z;imbdp=z;mask=false(siz);mask_shrinked=false(siz);
corr_im=z;
mask_corr_im=false(siz);
% r1 = 10; r2 = 150; 
r1 = 50; r2 = 200; 
corr_threshold=1; 
% x = NaN(Nframe,1e3); 
% y = NaN(Nframe,1e3); 
for ii = 1:Nframe
    corr_im(:,:,ii) = imresize(r.data_dic_save.displacements(ii).plot_corrcoef_dic,[siz(1),siz(2)]); 
end
mask_corr_ref = corr_im > corr_threshold; 
%
for ii = 1:Nframe
    u = imresize(r.data_dic_save.displacements(ii).plot_u_ref_formatted,[siz(1),siz(2)]);
    v = imresize(r.data_dic_save.displacements(ii).plot_v_ref_formatted,[siz(1),siz(2)]);
    [y0,x0] = find(mask_corr_ref(:,:,ii) == true);
    idx0 = sub2ind([siz(1),siz(2)],y0,x0);
    y=floor(y0+v(idx0)); x=floor(x0+u(idx0));
    idx = sub2ind([siz(1),siz(2)],y,x);
    mask_corr_im_ii=false([siz(1),siz(2)]);
    mask_corr_im_ii(idx)=true; 
    mask_corr_im(:,:,ii)=mask_corr_im_ii; 
end
%player(mask_corr_im);
%
for ii = 1:Nframe
    im_ii = satur(r.current_save(ii).gs*255,'level',70);
    mask(:,:,ii) = r.current_save(ii).roi.mask;
    mask_shrinked(:,:,ii) = shrink_mask(mask(:,:,ii),0.9);
    if ii == 1
        imbdp_ii = bandpassfft(im_ii,r1,r2,1);
    else
        imbdp_ii = bandpassfft(im_ii,r1,r2);
    end
    im_ii(~mask(:,:,ii)) = 0;
    imbdp_ii(~mask(:,:,ii)) = 0;
    im(:,:,ii) = reshape(im_ii,siz(1),siz(2));
    imbdp(:,:,ii) = reshape(imbdp_ii,siz(1),siz(2));
    if ii == 1
        y = prctile(imbdp_ii(mask(:,:,ii)),[1 99]);
    end
end
imbdp_norm10 = im2uint8((imbdp-y(1))/(y(2)-y(1)));

%%
player(imbdp_norm10);

%%
imbdp_norm_modif = imbdp_norm10;
% mask_corr = corr_im > 0.75; 
mask_bright = imbdp_norm_modif > 240; 
imbdp_norm_modif(mask_corr_im&mask_bright&mask_shrinked) = 50;
player(imbdp_norm_modif)

%%
player([mask_corr_im,mask_bright,mask_shrinked,mask_corr_im&mask_bright&mask_shrinked]);

%%

%% second filtering over modified image
imbdp2=z;
r1 = 10; 
r2 = 200; 
for ii = 1:Nframe
    im_ii = imbdp_norm_modif(:,:,ii);
    imbdp2_ii = bandpassfft(im_ii,r1,r2);
    imbdp2_ii(~mask(:,:,ii)) = 0;
    imbdp2(:,:,ii) = reshape(imbdp2_ii,siz(1),siz(2));
    if ii == 1
        y = prctile(imbdp_ii(mask(:,:,ii)),[1 99]);
    end
end
imbdp2_norm = im2uint8((imbdp2-y(1))/(y(2)-y(1)));
%%
player(imbdp2_norm)

%%
player([imbdp_norm1,imbdp_norm2,imbdp_norm3,imbdp_norm4,imbdp_norm5;...
    imbdp_norm6,imbdp_norm7,imbdp_norm8,imbdp_norm9,imbdp_norm10]);

%%
player(imbinarize(imbdp_norm))
%%
ims = im(:,:,1);
ims = ims(cursor_info(1).Position(2):cursor_info(2).Position(2),...
    cursor_info(1).Position(1):cursor_info(2).Position(1)); 

imfft = fft2(ims);%,2^nextpow2(1000),2^nextpow2(1000));

figure; 
ax = subplot(2,2,1);  
imshow(ims,[]);
ax = subplot(2,2,[2 4]);  
imhist(uint8(ims));  xlim([0 70]); ylim(ax,[0 200]); 
ax = subplot(2,2,3); hold(ax,'on');  
imshow(abs(fftshift(imfft)),[]); colormap(ax,'jet'); colorbar(gca);
caxis([0 1e3]); 
axis on 
%%
in = im(:,:,1); 
% Parameters
rmin = 10;
rmax = 500;
rthresh = rmax + 20;

% Scale input image between 0 and 1
I = imscale(in);

% Compute 2D DFT of input image
F = fft2(ifftshift(I));

% Extract magnitude
M = fftshift(F); % Center DFT
M = abs(M);      % Get the magnitude
M = log(M+1);    % Use log, for perceptual scaling, and +1 since log(0) is undefined
M = 100*M;       % Arbitrary scaling factor to improve image display
M = imscale(M);  % Scale the image between 0 and 1

% Remove DC and HF and take relevant image subset
siz = size(M);
c = ceil(siz/2);
xgv = (1:siz(2)) - c(2);
ygv = (1:siz(1)) - c(1);
[x,y] = meshgrid(xgv,ygv);
radius = sqrt(x.^2+y.^2);
M(radius < rmin) = 0;
M(radius > rmax) = 0;
M = M(c(1)-rthresh:c(1)+rthresh,c(2)-rthresh:c(2)+rthresh);

% Rescale image between 0 and 1
M = imscale(M);

figure; 
subplot(1,2,1); imshow(in,[]); 
subplot(1,2,2); imshow(M); colormap(gca,'jet'); colorbar(gca);