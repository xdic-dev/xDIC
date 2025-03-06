
%%% Effect of radius on MATCHING process 
% WITHOUT Reflection removal technics 
%% Load 
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 
subject = "S01"; 
material = "coating";
trial = "001";
phase = "loading"; 
ncorr1 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr12_40.mat'));
ncorr2 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr12_50.mat')); 
ncorr3 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr12_60.mat'));
ncorr4 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr12_70.mat'));

u = cat(3,ncorr1.data_dic_save.displacements(1).plot_u_ref_formatted,...
    ncorr2.data_dic_save.displacements(1).plot_u_ref_formatted,...
    ncorr3.data_dic_save.displacements(1).plot_u_ref_formatted,...
    ncorr4.data_dic_save.displacements(1).plot_u_ref_formatted); 

v = cat(3,ncorr1.data_dic_save.displacements(1).plot_v_ref_formatted,...
    ncorr2.data_dic_save.displacements(1).plot_v_ref_formatted,...
    ncorr3.data_dic_save.displacements(1).plot_v_ref_formatted,...
    ncorr4.data_dic_save.displacements(1).plot_v_ref_formatted); 

%%
newfig('displacement'); 
subplot(2,1,1); montage(u,'Size',[1 4],'BorderSize',0,'DisplayRange',[])
colormap(gca,'jet'); colorbar 
subplot(2,1,2); montage(v,'Size',[1 4],'BorderSize',0,'DisplayRange',[-75,-55])
colormap(gca,'jet'); colorbar
newfig('diff'); 
subplot(2,1,1); montage(u-u(:,:,4),'Size',[1 4],'BorderSize',0,'DisplayRange',[-2 2])
colormap(gca,'jet'); colorbar
subplot(2,1,2); montage(v-v(:,:,4),'Size',[1 4],'BorderSize',0,'DisplayRange',[-2 2])
colormap(gca,'jet'); colorbar

% We can see that the larger the radius, the less sensitive the matching is
% to reflections artefacts. It is what we expect. 
% The problem is that we see some smoothing of the field meaning that subsets 
% are no longer appropriatly approximated by linear first order
% transformations.

%% Effect of filtering on tracking process
%%%
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 
subject = "S01"; 
material = "coating";
trial = "001";
phase = "loading"; 
ncorr1 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'trial_with_filteredfield','radius40','ncorr1.mat'));
% ncorr2 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr1_raw_35.mat'));
% ncorr3 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr3.mat'));
% ncorr4 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr1_reflremov_30.mat'));
%%
% ncorr5 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'ncorr2.mat'));

%%
baseDataPath = "E:\transfert_SMP"; 
[cam_raw_1_all,cam_raw_2_all] = import_vid(baseDataPath,...
                   'subject',subject,...
                   'material',material,...
                   'trial',trial,...
                   'stereopair',1,...
                   'phase',phase,...
                   'framejump',1);
frameidx = load(fullfile(baseResultPath,'analysis',subject,'frameidx.mat')); 
frameidx = frameidx.frameidx; 
%% filtering images : use all images
% Those images were used for the generation of ncorr results 
limit_grayscale_level = 100; 
cam_bin_1_all = image_filtering(cam_raw_1_all,limit_grayscale_level); 
cam_reflremov_1_all = filteringtest(cam_raw_1_all);

jump = 2; 
%Keep only the targetted images 
startframe = frameidx{2}.loading.start(str2double(trial));
endframe = frameidx{2}.loading.end(str2double(trial));
cam_raw_1 = cam_raw_1_all(:,:,startframe:jump:endframe); 
cam_bin_1 = cam_bin_1_all(:,:,startframe:jump:endframe); 
cam_reflremov_1 = cam_reflremov_1_all(:,:,startframe:jump:endframe); 


cam_bin_2_all = image_filtering(cam_raw_2_all,limit_grayscale_level); 
cam_reflremov_2_all = filteringtest(cam_raw_2_all);

%Keep only the targetted images 
startframe = frameidx{2}.loading.start(str2double(trial));
endframe = frameidx{2}.loading.end(str2double(trial));
cam_raw_2 = cam_raw_2_all(:,:,startframe:jump:endframe); 
cam_bin_2 = cam_bin_2_all(:,:,startframe:jump:endframe); 
cam_reflremov_2 = cam_reflremov_2_all(:,:,startframe:jump:endframe); 
%%
player(cam_raw);player(cam_bin);player(cam_reflremov); 
%%
siz = size(cam_raw); 
hf = newfig; 
ax = subplot_ax(2,2); 
set(hf,'CurrentAxes',ax(1)); 
imshow(satur(cam_raw(:,:,end),'level',75),[]);  colormap(gca,'jet'); colorbar(gca);  
hold(ax,'on'); 
set(hf,'CurrentAxes',ax(2)); 
imshow(cam_reflremov(:,:,end),[]);  colormap(gca,'jet'); colorbar(gca);  
set(hf,'CurrentAxes',ax(3)); 
imshow(imresize(u1(:,:,end),[siz(1),siz(2)]),[]); colormap(gca,'jet'); colorbar(gca);  
set(hf,'CurrentAxes',ax(4)); 
imshow(imresize(u4(:,:,end),[siz(1),siz(2)]),[]); colormap(gca,'jet'); colorbar(gca);  
%%
player(satur(cam_raw,'level',75));colormap(gca,'jet'); colorbar(gca); 
%%

%% jet video masked by ncorr results 
siz = size(cam_raw_1); 
cam_mask = cam_raw_1; 
for ii = 1:siz(3)
    if siz(3)~= length(ncorr1.current_save)
        if ii == 1
            mask_ii = ncorr3.current_save(1).roi.mask; 
        else
            mask_ii = ncorr3.current_save(ii-1).roi.mask; 
        end
        im = cam_mask(:,:,ii); 
%         im(~mask_ii) = 0; 
        cam_mask(:,:,ii) = im; 
    else
        mask_ii = ncorr1.current_save(ii).roi.mask; 
        im = cam_mask(:,:,ii); 
        im(~mask_ii) = 0; 
        cam_mask(:,:,ii) = im; 
    end
end
%%
cam_mask4 = cam_mask; 
%%
player(satur(cam_mask,'level',90));
colormap(gca,'jet');colorbar(gca); 

%%
player([satur([cam_mask4],'level',90)],'movie',1);
colormap(gca,'jet');colorbar(gca); 
%%
imarray = cam_raw4; 
siz = size(imarray); 
% Filter image array with bandpassfft
r1 = 20; r2 = 200;
imbdp = zeros(siz);
y_all = zeros(1,2); 
for ii = 1:siz(3)
    imbdp(:,:,ii) = bandpassfft(imarray(:,:,ii),r1,r2);
    imbdp_ii = imbdp(:,:,ii);
    mask_ii = ncorr1.current_save(ii).roi.mask;
    % Remove extreme values from images
    y = prctile(imbdp_ii(mask_ii),[1 99]);
    y_all(1) = y_all(1)+y(1); 
    y_all(2) = y_all(2)+y(2); 
end
y_all = y_all/siz(3);
imbdp_norm = im2uint8((imbdp-y_all(1))/(y_all(2)-y_all(1)));

%%
player(imbdp_norm); 
%%
player(satur(imbdp,'level',65)); 

%%
player(cam_mask); 
%% sum histogram of masked video 
limit_x = 80; 
c = []; bin = []; 
% for ii = 1:siz(3)
    [counts,binLocations] = imhist(cam_mask); 
%     c = c+counts; 
%     bin = bin+binLocations; 
% end
binLocations = binLocations(2:end); 
counts = counts(2:end); 
fig = newfig; 
ax = subplot_ax(1,1); 
bar(ax(1),binLocations,counts); 
xlim(ax(1),[0,limit_x]); 
%%
Nframe = length(ncorr1.data_dic_save.displacements);
siz = size(ncorr1.data_dic_save.displacements(end).plot_u_ref_formatted);
z = zeros(siz(1),siz(2),Nframe); 
u1 = z; u2 = z; u3 = z;  u4 = z;
v1 = z; v2 = z; v3 = z;  v4 = z;
mask = z; 
for ii = 1:Nframe
    u1(:,:,ii) = ncorr1.data_dic_save.displacements(ii).plot_u_ref_formatted; 
%     u2(:,:,ii) = ncorr2.data_dic_save.displacements(ii).plot_u_ref_formatted; 
%     u3(:,:,ii) = ncorr3.data_dic_save.displacements(ii).plot_u_ref_formatted; 
%     u4(:,:,ii) = ncorr4.data_dic_save.displacements(ii).plot_u_ref_formatted; 
    
    v1(:,:,ii) = ncorr1.data_dic_save.displacements(ii).plot_v_ref_formatted; 
%     v2(:,:,ii) = ncorr2.data_dic_save.displacements(ii).plot_v_ref_formatted; 
%     v3(:,:,ii) = ncorr3.data_dic_save.displacements(ii).plot_v_ref_formatted; 
%     v4(:,:,ii) = ncorr4.data_dic_save.displacements(ii).plot_v_ref_formatted; 
    
%     mask(:,:,ii) = imresize(ncorr4.current_save(ii).roi.mask,[siz(1),siz(2)]);
    mask(:,:,ii) = ncorr1.data_dic_save.displacements(ii).roi_ref_formatted.mask;
end

%% visu diplacement
player([u1,u2,u3,u4]); colormap(gca,'jet'); colorbar 
player([v1,v2,v3,v4]); colormap(gca,'jet'); colorbar 
%% visu differences from the refl1remov
player([u1-u3,u2-u3,u3-u4],'imposedImLimit',[-1 1]); colormap(gca,'jet'); colorbar(gca);  
player([v1-v3,v2-v3,v3-v4],'imposedImLimit',[-1 1]); colormap(gca,'jet'); colorbar(gca);  
%% spatio-temporal low pass filter to displacement
siz = size(u1);
mask_corr_coef = zeros(siz);
level_corr_coef = 1.0; 
for ii = 1:siz(3)
    mask_corr_coef(:,:,ii)=ncorr1.data_dic_save.displacements(ii).plot_corrcoef_dic; 
    M(ii) = max(max(mask_corr_coef(:,:,ii)));
end
mask_corr_coef = logical(mask_corr_coef>level_corr_coef); 
u1modif = u1; v1modif = v1; 
% u1modif(mask_corr_coef) = NaN; 
% v1modif(mask_corr_coef) = NaN; 

% fill the NaN value
% u1modif = reshape(permute(u1modif,[3 1 2]),[siz(3),siz(1)*siz(2)]);
% v1modif = reshape(permute(u1modif,[3 1 2]),[siz(3),siz(1)*siz(2)]);
% v1modif = fillmissing(v1modif,'movmedian',10);  
% v1modif = fillmissing(v1modif,'movmedian',10);  
%DIC3DpairResults.corrComb{ii}=max([CorCoeff{ii} CorCoeff{ii+nImages}],[],2);

sigma = [1 1 2]; 
u1f = imgaussfilt3(u1modif, sigma,'FilterDomain','spatial','padding','replicate');
% u2f = imgaussfilt3(u2, sigma,'FilterDomain','spatial','padding','replicate');
% u3f = imgaussfilt3(u3, sigma,'FilterDomain','spatial','padding','replicate');
% u4f = imgaussfilt3(u4, sigma,'FilterDomain','spatial','padding','replicate');

v1f = imgaussfilt3(v1modif, sigma,'FilterDomain','spatial','padding','replicate');
% v2f = imgaussfilt3(v2, sigma,'FilterDomain','spatial','padding','replicate');
% v3f = imgaussfilt3(v3, sigma,'FilterDomain','spatial','padding','replicate');
% v4f = imgaussfilt3(v4, sigma,'FilterDomain','spatial','padding','replicate');

threshold=10;
% mask_diff=~(abs(v1-v1f)>threshold); 
mask_in = u1(:,:,end)~=0; 
[mask_out] = shrink_mask(mask_in,0.9);
mask_in=repmat(mask_in,[1 1 Nframe]);
mask_out=repmat(mask_out,[1 1 Nframe]);
player([mask_in,mask_out,mask_corr_coef&mask_out]);
% u1modif(mask_corr_coef&mask_out)=u1f(mask_corr_coef&mask_out);
% v1modif(mask_corr_coef&mask_out)=v1f(mask_corr_coef&mask_out);
player([v1,mask_corr_coef.*max(v1,[],'all'),v1modif,v1f]);colormap(gca,'jet'); colorbar(gca);  

%%

u1modif(mask_out)=u1(mask_corr_coef&mask_out);
v1modif(mask_out)=v1f(mask_corr_coef&mask_out);
player(u1modif);colormap(gca,'jet'); colorbar(gca);  
%%
% player([u1-u1f,v1-v1f,mask_corr_coef]); colormap(gca,'jet'); colorbar(gca);  
player([u1,u1modif]); colormap(gca,'jet'); colorbar(gca);  
% player([v1,v1modif]); colormap(gca,'jet'); colorbar(gca);  
%%
player([v1,v1modif]); colormap(gca,'jet'); colorbar(gca);  

%player(v1f); colormap(gca,'jet'); colorbar(gca);  
%%
level = 0.5; 
masku1 = abs(u1f - u1)>level; 
masku2 = abs(u2f - u2)>level; 
masku3 = abs(u3f - u3)>level; 
masku4 = abs(u4f - u4)>level; 

maskv1 = abs(v1f - v1)>level; 
maskv2 = abs(v2f - v2)>level; 
maskv3 = abs(v3f - v3)>level; 
maskv4 = abs(v4f - v4)>level; 
%%
player([masku1],'movie',1,'fps',10); %,masku2,masku3,masku4]); 
%%
player([u1,u2,u3,u4;u1f,u2f,u3f,u4f]); colormap(gca,'jet'); colorbar(gca);  

%% visu differences from the refl1remov
player([abs(u1-u1f)],'imposedImLimit',[0 1],'movie',1,'fps',10); colormap(gca,'jet'); colorbar(gca);  
% player([u1-u1f,u2-u2f,u3-u3f,u4-u4f],'imposedImLimit',[-1 1]); colormap(gca,'jet'); colorbar(gca);  
% player([v1-v1f,v2-v2f,v3-v3f,v4-v4f],'imposedImLimit',[-1 1]); colormap(gca,'jet'); colorbar(gca);  

%%
u4_modif = u4; 
v4_modif = v4; 
for ii = 3:Nframe-2 
   tempregionmask = mask(:,:,ii); 
   [tempregionmask] = shrink_mask(tempregionmask, 0.6);
   
   tempu = u4_modif(:,:,ii); 
   tempuf = u4f(:,:,ii); 
   tempumask = masku4(:,:,ii); 
   tempumask = tempregionmask & tempumask; 
   
   tempv = v4_modif(:,:,ii); 
   tempvf = v4f(:,:,ii); 
   tempvmask = maskv4(:,:,ii); 
   tempvmask = tempregionmask & tempvmask; 
   
   tempu(tempumask) = tempuf(tempumask);%use filter data only at mask location
   tempv(tempvmask) = tempvf(tempvmask);
   
   u4_modif(:,:,ii) = tempu; 
   v4_modif(:,:,ii) = tempv; 
    
end
%%
figure; imshow(tempumask); 
%%
player(abs(u4-u4_modif));colormap(gca,'jet'); colorbar(gca);
%%
player([u4,u4_modif;masku4*15,masku4*15]);colormap(gca,'jet'); colorbar(gca);

%% 
%%% Filtering displacement field from saturated images (number 1): 
% First find the subset that were badly tracked. (already done in a
% previous section).

%%
u1_modif = u1; 
v1_modif = v1; 
for ii = 3:Nframe-2 
   tempregionmask = mask(:,:,ii); %mask of the image
   [tempregionmask] = shrink_mask(tempregionmask, 0.6);
   
   tempu = u1_modif(:,:,ii); 
   tempuf = u1f(:,:,ii); 
   tempumask = masku4(:,:,ii); 
   tempumask = tempregionmask & tempumask; 
   
   tempv = v1_modif(:,:,ii); 
   tempvf = v1f(:,:,ii); 
   tempvmask = maskv1(:,:,ii); 
   tempvmask = tempregionmask & tempvmask; 
   
   tempu(tempumask) = tempuf(tempumask);%use filter data only at mask location
   tempv(tempvmask) = tempvf(tempvmask);
   
   u1_modif(:,:,ii) = tempu; 
   v1_modif(:,:,ii) = tempv;  
end


%%
field = u1; 
siz = size(field); 
fieldperm = permute(u1,[3 1 2]); 
fieldcol = reshape(fieldperm,siz(3),siz(1)*siz(2),1);

idx = ~any(fieldcol==0);
masknan = ~any(isnan(fieldcol)); 
fieldm = fieldcol(:,idx&masknan);
data = fieldm;
data = [data(1,:)-data(1,:);diff(data)];
npad = 15; 
freqFilt = 4; %number of image over the filter is performed.
freqAcq = siz(3);
x = [repmat(data(1,:),[npad,1]);...
    data(:,:);...
    repmat(data(end,:),[npad,1])];
[B,A] = butter(4,freqFilt/(freqAcq/2)); 
data_f = filtfilt(B,A,x);
data_f = data_f(npad:end-npad-1,:);
data_f = fieldm(1,:)+cumsum(data_f);
data_f_clean = zeros(siz(3),length(idx)); 
data_f_clean(:,idx&masknan) = data_f;

fieldperm_f = reshape(data_f_clean,siz(3),siz(1),siz(2));
field_f = permute(fieldperm_f,[2 3 1]); 

%player(field_f-field);colormap(gca,'jet'); colorbar(gca);
figure; 
imshow(data_f(:,400:500)); colormap(gca,'jet'); colorbar(gca);
%%
player([field,field_f,u1f],'imposedImLimit',[-20,25],'movie',1,'fps',10);colormap(gca,'jet'); colorbar(gca);
%%
siz = size(u1f); 
u1fcol = permute(u1f,[3 1 2]); 
u1fcol = reshape(u1fcol,siz(3),siz(1)*siz(2),1);
a = u1fcol(:,idx); 

newfig('field');
ax = subplot_ax(2,3); hold(ax,'on');
plot(ax(1),1:siz(3),fieldm,'.-');
plot(ax(2),1:siz(3),data_f,'.-');
plot(ax(3),1:siz(3),a,'.-');
plot(ax(4),1:siz(3),[fieldm(1,:)-fieldm(1,:);diff(fieldm)],'.-'); ylim(ax(4),[-5,5]);
plot(ax(5),1:siz(3),[data_f(1,:)-data_f(1,:);diff(data_f)],'.-'); ylim(ax(5),[-5,5]);
plot(ax(6),1:siz(3),[a(1,:)-a(1,:);diff(a)],'.-'); ylim(ax(6),[-5,5]);
linkaxes(ax,'x');
% linkaxes([ax(1),ax(2)],'y'); 
% linkaxes([ax(3),ax(4)],'y'); 
%%
data = fieldm(:,:);
data = [data(1,:)-data(1,:);diff(data)];
A = data; 
% A = [60 59 49 49 58 100 61 57 48 58];
[B,TFrm] = rmoutliers(A, 'percentiles', [5 95]);
figure; plot(B)
% figure; 
% plot(A,'r')
% hold on
% plot(find(~TFrm),B,"bo-")
% % yline([L U C],":",["Lower Threshold","Upper Threshold","Center Value"])
% legend("Original Data","Cleaned Data")

%%
im = satur(cam_raw_1,'level',70);
size(im) 
player(im)
%% Images reconstruction from ncorr results
% We see that overall the grayvalues are conserved.
siz = size(ncorr1.current_save(1).gs); 
Nframe = length(ncorr1.current_save); 
siz = [siz,Nframe]; 
z = zeros(siz); 
im_ref = z; im_cur = z; 
[y0,x0] = find(imresize(u1(:,:,1),[siz(1),siz(2)])); 
[idx] = find(imresize(u1(:,:,1),[siz(1),siz(2)])); %mask from non zeros deplacement 

mask_disp_ref = u1(:,:,1)~=0; 
mask_disp_ref_red = shrink_mask(mask_disp_ref,0.9);
mask_border = mask_disp_ref&~mask_disp_ref_red;

for ii = 1:Nframe
    %initialisation
    im_ref_ii = im_ref(:,:,ii);
    im_cur_ii = im_cur(:,:,ii);
    %real mask use in ncorr
    mask_im = ncorr1.current_save(ii).roi.mask;
    %remove border of displacement field filtered
    u1_ii = u1(:,:,ii); v1_ii = v1(:,:,ii); %unfiltered displacement field 
    u_ii = u1f(:,:,ii); v_ii = v1f(:,:,ii); 
    u_ii(mask_border) = u1_ii(mask_border); u_ii=reshape(u_ii,size(u1_ii)); 
    v_ii(mask_border) = v1_ii(mask_border); v_ii=reshape(v_ii,size(u1_ii)); 
    %interpolate between untracked pixels (subset untracked) 
    u_ii = imresize(u_ii,[siz(1),siz(2)]); 
    v_ii = imresize(v_ii,[siz(1),siz(2)]);  
    %pixel displacement 
    y_ii = y0+v_ii(idx);
    x_ii = x0+u_ii(idx);
    %indices of neighboring pixels 
    kk = 0;
    z = zeros(length(x_ii),4); 
    dist_ii=z; ind_ii=z; gs_ii=z; 
    y_kk = [floor(y_ii),ceil(y_ii)];
    x_kk = [floor(x_ii),ceil(x_ii)];
    for oo = 1:2
        for pp = 1:2
            kk = kk+1; 
            dist_ii(:,kk) = sqrt((y_ii-y_kk(:,oo)).^2+(x_ii-x_kk(:,pp)).^2);
            ind_ii(:,kk) = sub2ind([siz(1),siz(2)],y_kk(:,oo),x_kk(:,pp));
            gs_ii(:,kk) = ncorr1.current_save(ii).gs(ind_ii(:,kk));
        end
    end
    if ii == 1
        gs_averaged = ncorr1.current_save(ii).gs(sub2ind([siz(1),siz(2)],round(y0),round(x0)));
    else    
        w_ii = (1./dist_ii); %weight 
        w_norm_ii = w_ii./sum(w_ii,2); %weight normalized
        gs_averaged = sum(w_norm_ii.*gs_ii,2); %average over neighborhood 
    end
    
    idx_ii = sub2ind([siz(1),siz(2)],round(y_ii),round(x_ii)); 
    im_ref_ii(idx) = gs_averaged;
    im_cur_ii(idx_ii) = gs_averaged;
    
    %output in MxNxNframe
    im_ref(:,:,ii) = reshape(im_ref_ii,[siz(1),siz(2)]); 
    im_cur(:,:,ii) = reshape(im_cur_ii,[siz(1),siz(2)]); 
end
% sigma = [0.01 0.01 3]; 
% im_ref = imgaussfilt3(im_ref, sigma,'FilterDomain','spatial','padding','replicate');
% im_cur = imgaussfilt3(im_cur, sigma,'FilterDomain','spatial','padding','replicate');
%%
player([im_ref],'movie',1,'fps',10);
%%
player(satur([im_cur],'level',75/255),'movie',1,'fps',10);
%%
siz = size(ncorr1.current_save(ii).gs); 
Nframe = length(ncorr1.current_save);
siz = [siz,Nframe];
im = zeros(siz); 
for ii = 1:Nframe
  im_ii = ncorr1.current_save(ii).gs; 
  im_ii(~ncorr1.current_save(ii).roi.mask) = 0; 
  im(:,:,ii) = im_ii; 
end
player(satur(im,'level',75/255),'movie',1,'fps',10)
%%
mask = im_ref(:,:,1) == 0; 
Nvec = length(find(~mask));
im_vec = zeros(Nvec,Nframe); 
for ii =  1:Nframe
    im_ii = im_ref(:,:,ii); 
    im_ii(mask) = []; 
    im_vec(:,ii) = im_ii; 
end
%%
figure; 
plot(1:Nframe,im_vec(1:1000,:)');
ylim([0 1]); 
%%
mask3D = repmat(ncorr1.reference_save.roi.mask,[1 1 Nframe]);
imvec = im(mask3D);
imvec = reshape(imvec,[],length(imvec)/Nframe);
%%
imvec = reshape(permute(im,[3 1 2]),[Nframe siz(1)*siz(2)]);
mask = imvec == 0; 
imvec(mask) = NaN; 
imvec = reshape(imvec,Nframe,[]);
%%
mask_nan = isnan(imvec);
%%
figure; 
plot(1:Nframe,imvec)
%% mean image
figure; imshow(mean(im,3,'omitnan'));
%%
im = satur(cam_raw_1,'level',75); 
%% correlation
im_level_min = 255*ones(siz(1),siz(2)); 
for ii = 1:Nframe
    im_ii = im(:,:,ii); 
    R(ii) = corr2(im(:,:,1),im_ii);
    im_level_min(im_level_min>im_ii) = im_ii(im_level_min>im_ii); 
    im_level_min = reshape(im_level_min,[siz(1),siz(2)]);
end
figure; 
plot(1:Nframe,R,'.');
ylim([0.8 1]); 
ylabel('Correlation');
xlabel('frame #'); 
%%
figure; 
imshow(medfilt2(im(:,:,1)));
%%
level = 75; 
mask = im>level; 
im_modif = im;
player(satur(im,'level',level)); colormap(gca,'jet'); colorbar(gca);

%% mask with only level intensity criterium and graylevel intensity replaced 
% by gray scale found when 
level = 70;
mask = imarray>level; 
for ii = 1:Nframe
   tempregionmask = mask(:,:,ii); 
   [tempregionmask] = shrink_mask(tempregionmask, 0.6);
   
   tempu = u4_modif(:,:,ii); 
   tempuf = u4f(:,:,ii); 
   tempumask = mask(:,:,ii); 
   tempumask = tempregionmask & tempumask; 
   tempu(tempumask) = tempuf(tempumask);%use filter data only at mask location
   im_modif(:,:,ii) = tempu; 
    
end
%% addition useful functions 
function [mask_out,B] = shrink_mask(mask, factor)
siz = size(mask); 
B = bwboundaries(imgradient(mask),'holes');
A = B{1}; 
xi = A(:,1);
yi = A(:,2);
j = boundary(xi,yi); 
xi = xi(j);
yi = yi(j);
p = polyshape(xi,yi);
[cx, cy] = centroid(p); %find the centroid
m = scale(p,factor,[cx cy]);
xi = m.Vertices(:,1);
yi = m.Vertices(:,2);
mask_out = logical(poly2mask(yi,xi,siz(1),siz(2)));
end