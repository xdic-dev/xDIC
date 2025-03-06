

%%% ncorr JAI analysis 
baseDataPath = "F:\transfert_SMP"; 
baseResultPath = "E:\data_dic\fingertip3drecontruction"; 
subject = "S09"; % S03 
material = "coating";
trial = "013"; %025
phase = "loading"; 

% Visu trial
frameStart = 1; 
frameEnd = 150; 
framejump = 1; 
trialname = subject+"_"+material+"_speckles_"+trial; 
%ftir
file = fullfile(baseDataPath,'rawdata',subject,'speckles',material,'vid',sprintf('%s*_cam_0_0.mp4', trialname)); 
S = dir(file);
im_fitr_raw = readvid(fullfile(S.folder, S.name),(frameEnd-frameStart+1)/framejump,frameStart-1,framejump);
%speckles
file = fullfile(baseDataPath,'rawdata',subject,'speckles',material,'vid',sprintf('%s*_cam_0_1.mp4', trialname));
S = dir(file);
im_speckles_raw = readvid(fullfile(S.folder, S.name),(frameEnd-frameStart+1)/framejump,frameStart-1,framejump);

%% visu 
player([im_fitr_raw,im_speckles_raw],'movie',0);
%%
im_test = -((im_speckles_bp(:,:,span))-(im_fitr_bp(:,:,span)));
% im_test_1 = im_test(:,:,end/2); 
% player(satur(im_test,'method','low','level',240));
% player(im_test,'contour',masks_saved.caxy_fitr);%,im_test]); 
  
y = prctile(im_test(:,:,end),[0.5 99.5],'all');
im_test = im2double((im_test-y(1))/(y(2)-y(1)));

%%
player(im_speckles_bp)
%%
player(im_fitr_bp)
%%

%%
player(im_test,'contour',masks_saved.caxy_fitr);
  
%%
player([im_speckles_bp,im_fitr_bp,im_test],'contour',masks_saved.caxy_fitr+[2*size(im_test,2),0]); 
  
%% overall param
id_start = 8;%13; 
id_end = 100+id_start-1;%113-1;
jump=1; 
im_fitr = im_fitr_raw(:,:,id_start:jump:id_end);
im_speckles = im_speckles_raw(:,:,id_start:jump:id_end);

% Remove extreme values from images
y = prctile(im_fitr(:,:,1:20:end),[0.5 99.5],'all');
im_fitr_m = im2double(im_fitr-y(1)/(y(2)-y(1)));
%%
y = prctile(im_speckles(:,:,1:20:end),[0.5 99.5],'all');
im_speckles_m = im2double(im_speckles-y(1)/(y(2)-y(1)));

% Attempt to remove speckles from fitr images
BW = false(size(im_fitr_m)); 
for ii = 1:size(im_fitr,3)
    [BW(:,:,ii)] = imbinarize(im_speckles_m(:,:,ii),'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
end
im_fitr_m(~BW) = 0.4; %mean(imf,'all'); 
%%
im_modif = im_fitr_m(:,:,frame_span); 
%%
player([im_fitr_m(:,:,frame_span),im_modif],'contour',[masks_saved.caxy_fitr(:,:,frame_span);masks_saved.caxy_fitr(:,:,frame_span)+[size(im_fitr(:,:,1),2),0]],'movie',1); 
%%
player(im_speckles_m(:,:,frame_span),'movie',1)
%% Run benoit's method 
%% bandpass filter 
% hpf=padarray(1,[5 5])-fspecial('gaussian',11,2);
% lpf=fspecial('gaussian',11,2);
% im_fitr_bp=imfilter(imfilter(im2single(im_fitr),hpf,'replicate'),lpf,'replicate');
% im_fitr_bp = satur(satur(im_fitr,'method','high','level',85),'method','low','level',65); 
im_speckles_bp = zeros(size(im_fitr)); 
for ii = 1:size(im_fitr,3)
    im_speckles_bp(:,:,ii) = bandpassfft(im_speckles(:,:,ii),20,120);
%     im_fitr_bp(:,:,ii)=bandpassfft(satur(satur(im_fitr(:,:,ii),'method','high','level',100),'method','low','level',70),100,150);
%     if ii == size(im_fitr,3)
%         bandpassfft(im_fitr_bp(:,:,ii),20,400,1);
%     end
end
%%
% Remove extreme values from images
y = prctile(im_speckles(:,:,end),[0.5 99.5]);
im_speckles_bp = im2uint8((im_speckles-y(1))/(y(2)-y(1)));
%%
player([im_fitr_bp-(255-im_speckles_bp)]); 
%%
player(cat(3,diff(im_fitr(:,:,1:2),1,3),diff(im_fitr,1,3)))
%%
player([im_fitr(:,:,4:end),im_fitr_bp(:,:,4:end)],'contour',cat(1,masks_saved.caxy_fitr,masks_saved.caxy_fitr+[size(im_fitr,2),0]))
%%
player(cat(2,im_fitr_bp,im_fitr),'movie',0,'fps',10);

%% manual CA : manual draw the contour of the contact 
jumpCA = 20; 
idxCA = [5,7,9,13,18,24,30,40,70,100]; 
masks = manualCA(im_fitr(:,:,idxCA)); 

%% train and compute the contour of each image twice
[~,~,mdl] = AT_CA(im_fitr(:,:,5:end), masks, idxCA-4); 
[cabw_fitr,caxy_fitr] = AT_CA(im_fitr(:,:,5:end),[],[],mdl); 

%% save contour from benoit's method
save(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'mask_fitr_jai.mat'),'cabw_fitr','caxy_fitr');

%% load mask
masks_saved = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,'mask_fitr_spremoval_jai.mat')); 
sizMask = size(masks_saved.caxy_fitr);
masks_saved.caxy_fitr = cat(3,NaN(sizMask(1),sizMask(2),size(im_fitr_m,3)-sizMask(3)),masks_saved.caxy_fitr);
sizMask = size(masks_saved.cabw_fitr);
masks_saved.cabw_fitr = cat(3,NaN(sizMask(1),sizMask(2),size(im_fitr_m,3)-sizMask(3)),masks_saved.cabw_fitr);

%%
k=boundary(caxy_fitr(:,1,ii),caxy_fitr(:,2,ii),0);

%% visu contour
player([im_fitr_m,im_speckles_m],'contour',[masks_saved.caxy_fitr;masks_saved.caxy_fitr+[size(im_fitr(:,:,1),2),0]]); 

%% features tracking : find the best features to track and execute 
% the algorithm over all images 
tic; 
flxy_fitr = at_camera_OF_midpoint(...
    im_fitr(:,:,4:end),... 
    cabw(:,:,[1 round(end/2) end]),... %roi for optical flow
    [],... %subregion to track
    'spacing',30); % min distance in pxl between two features
toc; 
%%
im_speckles = im_speckles_raw; 
y = im2double(prctile(im_speckles(:,:,end),[1 99],'all'));
im_speckles_bp = (im2double(im_speckles)-y(1))/(y(2)-y(1));
im_speckles_bp = satur(satur(im_speckles_bp,'method','low','level',0),...
    'method','high','level',1);
player(im_speckles_bp)

%% Run dic analysis 
h = ncorr;
h.set_ref(im_speckles_bp(:,:,1)); 
h.set_cur(Myarray2cell(im_speckles_bp)); 
% h.set_roi_cur(masks_saved.cabw_fitr(:,:,end));

%% save dic analysis 
filename_dic = fullfile(baseResultPath,subject,material,trial,phase,'dic_jai.mat'); 
save_dicdata(h, filename_dic);

%% Load dic analysis 
% filename_dic = fullfile(baseResultPath,subject,material,trial,phase,'dic_1-4_jai.mat'); 
% ncorr1 = load(filename_dic);
% ncorrJAI = merge_ncorr_results(ncorr1,ncorr2); 
filename_dic = fullfile(baseResultPath,subject,material,trial,phase,'dic_jai.mat'); 
ncorrJAI = load(filename_dic);

%% Exploit dic tracking to extract strains 
sizIM = size(ncorrJAI.current_save(1).gs); 
Nframe = length(ncorrJAI.current_save); 
sizIM = [sizIM,Nframe]; 
fps = 50; 

sizDisp = size(ncorrJAI.data_dic_save.displacements(1).plot_u_ref_formatted); 
sizDisp = [sizDisp,Nframe]; 
z = zeros(sizDisp); 
u = z; v = z; mask = z; 
z = zeros(sizIM); 
im = z; cabw = imrotate(false(sizIM),0); 
caxy = zeros(1e3,2,Nframe); 
for ii = 1:Nframe 
    u(:,:,ii) = ncorrJAI.data_dic_save.displacements(ii).plot_u_ref_formatted; 
    v(:,:,ii) = ncorrJAI.data_dic_save.displacements(ii).plot_v_ref_formatted; 
    im(:,:,ii) = ncorrJAI.current_save(ii).gs;
    cabw(:,:,ii) = ncorrJAI.current_save(ii).roi.mask;
end

mask = u~=0; 
[y0,x0] = find(mask(:,:,1)~=NaN); 
Nfeat = length(reshape(mask(:,:,1),[],1));
flxy = NaN(Nfeat,2,Nframe); 

factorx = sizIM(2)/sizDisp(2); 
factory = sizIM(1)/sizDisp(1);

for ii = 1:Nframe
    u_ii = u(:,:,ii); v_ii = v(:,:,ii);% mask_ii = reshape(mask(:,:,ii),[],1);
    flxy(:,1,ii) = x0.*factorx + reshape(u_ii,[],1); 
    flxy(:,2,ii) = y0.*factory + reshape(v_ii,[],1); 
    mask(:,:,ii) = u_ii~=0; 
end
flxy(any(mask==0,3),:,:) = []; 
% flxyflip = flxy; 
% flxyflip(:,1,:) = flxy(:,2,:); 
% flxyflip(:,2,:) = flxy(:,1,:); 

% filter 
flxyfilt = filter_feature_position(flxy);

st = at_camera_ST(cabw,[sizIM(1),sizIM(2),sizIM(3)],flxyfilt,'smoothTime',true,'smoothSpace',true); 
st = adden(st); % add absolute value
st = cat(2,st,strain_pc(st(:,3+1,:),st(:,3+2,:),st(:,3+3,:))); % add principal components 
% st(:,[4,5,6,8,9,10],:) = st(:,[4,5,6,8,9,10],:)*fps*1e2; % to put in [%/s]
st(:,4:end,:)=st(:,4:end,:)*50*1e2; 
stcum = st; 
stcum(:,4:end,:) = cumsum(stcum(:,4:end,:),3,'omitnan')/50;
stcum(:,4:end,:) = stcum(:,4:end,:)-stcum(:,4:end,27);  
disp('done'); 
%%
[e,v] = strain_pc(st(:,3+1,:),st(:,3+2,:),st(:,3+3,:));
%%
frame_span = 1:99; 
%% exx 
player(im(:,:,frame_span),'contour',masks_saved.caxy_fitr(:,:,frame_span),'points',flxyfilt(:,:,frame_span),'strain',st(:,[1 2 4],frame_span),'vector',vector_displ_true(:,:,frame_span),'fps',50,'movie',0,'flipcmap',false,'imposedStLimit',0.005); 
%% eyy
player(im(:,:,frame_span),'contour',masks_saved.caxy_fitr(:,:,frame_span),'points',flxyfilt(:,:,frame_span),'strain',st(:,[1 2 5],frame_span),'vector',vector_displ_true(:,:,frame_span),'fps',50,'movie',0,'flipcmap',false,'imposedStLimit',0.005); 
%% e1
player(im(:,:,:),'points',flxyfilt(:,:,:),'strain',stcum(:,[1 2 9],:),'fps',50,'movie',0,'flipcmap',false,'imposedStLimit',10); 
%% e2
player(im(:,:,:),'points',flxyfilt(:,:,:),'strain',stcum(:,[1 2 10],:),'fps',50,'movie',0,'flipcmap',false,'imposedStLimit',10); 
%%
strain_show(st,'contour',masks_saved.caxy_fitr);

%% Quiver plot like Willemet 
ref_id = 1;  
cur_id = 20; 

jump_point = 10; 
Vfirst = [flxyfilt(1:jump_point:end,1,ref_id),flxyfilt(1:jump_point:end,2,ref_id)]; %position
Dnow = [flxyfilt(1:jump_point:end,1,cur_id)-flxyfilt(1:jump_point:end,1,ref_id),flxyfilt(1:jump_point:end,2,cur_id)-flxyfilt(1:jump_point:end,2,ref_id)]; %end point

title_fig = "quiver_displacement field"; 
newfig(title_fig); 
ax = subplot_ax(1,1); hold(ax,'on'); 
quiver(ax(1),Vfirst(:,1),Vfirst(:,2),Dnow(:,1),Dnow(:,2),0,'Color',.2*[1 1 1],'ShowArrowHead','on','AutoScale','off'); hold on;

%%
jump_point = 1; 
ref_id = 1; 
touch_id = 1; 
vector_displ = zeros(size(flxyfilt,1)/jump_point,4,Nframe); 
vector_displ_true = vector_displ; 
Vfirst = [flxyfilt(1:jump_point:end,1,ref_id),flxyfilt(1:jump_point:end,2,ref_id)]; 
% cur_id = 20; 
% Dlast = [flxyfilt(1:jump_point:end,1,cur_id),flxyfilt(1:jump_point:end,2,cur_id)]-Vfirst; %end point
% [~,idx] = min(abs(Dlast(:,2)));
siz = size(im); 
% figure; 
% hold(gca,'on'); 
% imshow(im(:,:,end)); 
% h = drawpoint();
% idx_x = find_closest_value(f,val)h.Position(1); 

for ii = 1:Nframe 
    cur_id = ii; 
    Dcenternow = [flxyfilt(idx,1,touch_id),flxyfilt(idx,2,touch_id)]-Vfirst(idx,:);  
    Dnow = [flxyfilt(1:jump_point:end,1,cur_id),flxyfilt(1:jump_point:end,2,cur_id)]-Vfirst; %end point
    vector_displ(:,:,ii) = [Vfirst,Dnow-Dcenternow];
    vector_displ_true(:,:,ii) = [Vfirst,Dnow]; 
end
% player([im,im],'points',flxyfilt,'morepoints',flxyfilt(idx,:,:),'vector',cat(1,vector_displ2,vector_displ+[siz(2),0,0,0]),'flipcmap',false,'imposedStLimit',0.01,'movie',1); 
% player([im*255,im_fitr(:,:,4:end)],'contour',cat(1,masks_saved.caxy_fitr,masks_saved.caxy_fitr+[size(im,2),0]),'points',flxyfilt,'morepoints',flxyfilt(idx,:,:),'vector',vector_displ_true,'flipcmap',false,'imposedStLimit',0.01,'movie',0); 

%% MAXIMUM CONTRACTION 
player(im,'contour',masks_saved.caxy_fitr,'strain',st(:,[1 2 10],:),'flipcmap',false,'imposedStLimit',0.005,'movie',0); 


%,'strain',cat(2,st(:,[1 2],:)cumsum(st(:,4,:),3))
% player(im,'points',flxyfilt,'morepoints',flxyfilt(idx,:,:),'vector',vector_displ2,'strain',st(:,[1 2 10],:),'flipcmap',false,'imposedStLimit',0.05); 

%%
frame_span = 1:30; 
player(im(:,:,frame_span),'contour',masks_saved.caxy_fitr,'points',flxyfilt(:,:,frame_span),...
    'vector',cat(2,st(:,[1 2],frame_span),v(:,[3 4],frame_span)*10),...
    'flipcmap',false,'imposedStLimit',0.050,'strain',cat(2,st(:,[1 2],frame_span),cumsum(st(:,10,frame_span),3,'omitnan')),'fps',10,'movie',1);%,...
%% 
player(im(1242-20:1242+20,:,:),'contour',masks_saved.caxy_fitr,'points',flxyfilt-[0,1242],'strain',st(:,[1 2 4],:)-[0,1242,0],'flipcmap',false,'imposedStLimit',0.005);

%% understand and visualize the motion of material points 
t = 15; %threshold 
title_fig = "Radial point trajectory without dephasing";
newfig(title_fig); 
ax = subplot_ax(1,2); hold(ax,'on'); 
scattersize = 25; 
time = linspace(2,Nframe,Nframe-1)/50; 
ref_id = 1; 
xcenter = 1313; ycenter =  1247; 
Vfirst = [st(:,1,ref_id),st(:,2,ref_id)]; 
resolution = 80; 
%
xval = flip(600:20:1800); yval = repelem(ycenter,length(xval)); 
xvalii = xcenter; yvalii = ycenter; 
idx_center = find((xvalii-t<=Vfirst(:,1)&Vfirst(:,1)<=xvalii+t)&...
    (yvalii-t<=Vfirst(:,2)&Vfirst(:,2)<=yvalii+t),1,'last');

% idx = zeros(length(xval),1);
idx = []; 
pp=0; 

for ii = 1:length(xval)
    xvalii = xval(ii);
    yvalii = yval(ii);
    idx_ii = find((xvalii-t<=Vfirst(:,1)&Vfirst(:,1)<=xvalii+t)&...
        (yvalii-t<=Vfirst(:,2)&Vfirst(:,2)<=yvalii+t),1,'last');
    ax_idx = 1; 
    m = false(Nframe-1,1); 
%     x = (squeeze(st(idx_ii,1,:))-squeeze(st(idx_center,1,:)))/resolution; 
    for jj = 5:Nframe-1
        m(jj) = masks_saved.cabw_fitr(round(squeeze(st(idx_ii,2,jj))),round(squeeze(st(idx_ii,1,jj))),jj);
    end
    x = (squeeze(st(idx_ii,1,:)))/resolution; 
    y = time;
    z = squeeze(st(idx_ii,4,:))*50*1e2;
    plot(ax(ax_idx),x,y,'k-'); 
    scatter(ax(ax_idx),x(~m),y(~m),...
        scattersize,...
        z(~m),'^','filled',...
        'LineWidth',1,...
        'MarkerEdgeColor',[1 1 1].*0); 
    scatter(ax(ax_idx),x(m),y(m),...
        scattersize,...
        z(m),'s','filled'); 
    

    
%     if mod(ii,10)==0 
%         ax_idx = 3; 
%         pp = pp+1; 
%         plot(ax(ax_idx),y,z)
%         idx(pp) = idx_ii; 
%     end
end
% ax.YTickLabel = num2str(flip(time'));
ax_idx = 1; 
colormap(ax(ax_idx),'jet'); 
cb=colorbar(ax(ax_idx)); 
xlabel(ax(ax_idx),'X pos [mm]');
ylabel(ax(ax_idx),'Time [sec]'); 
caxis(ax(ax_idx),[-10,10]); 
set(ax(ax_idx), 'YDir','reverse')
cb.Label.String = 'd/dt(e_{xx})';
set(ax(ax_idx),'FontSize',15); 
% ax_idx = 3; 
% legdata = (squeeze(st(idx,1,end))-squeeze(st(idx_center,1,end)))/resolution; 
% legend(ax(ax_idx),compose("%1.2f",legdata));
% set(ax(ax_idx), 'XDir','reverse')
%
yval = flip(700:20:1500); xval = repelem(xcenter,length(yval)); 
% idx = zeros(length(yval),1);
pp = 0; 
idx = []; 
for ii = 1:length(xval)
    ax_idx = 2; 
    xvalii = xval(ii);
    yvalii = yval(ii);
    idx_ii = find((xvalii-t<=Vfirst(:,1)&Vfirst(:,1)<=xvalii+t)&...
        (yvalii-t<=Vfirst(:,2)&Vfirst(:,2)<=yvalii+t),1,'last');
    for jj = 5:Nframe-1
        m(jj) = masks_saved.cabw_fitr(round(squeeze(st(idx_ii,2,jj))),round(squeeze(st(idx_ii,1,jj))),jj);
    end
    x = (squeeze(st(idx_ii,2,:))-ycenter)/resolution; 
    y = time; 
    z = squeeze(st(idx_ii,5,:))*50*1e2; 
    plot(ax(ax_idx),x,y,'k-'); 
    scatter(ax(ax_idx),x(~m),y(~m),...
        scattersize,...
        z(~m),'^','filled',...
        'LineWidth',1,...
        'MarkerEdgeColor',[1 1 1].*0); 
    scatter(ax(ax_idx),x(m),y(m),...
        scattersize,...
        z(m),'s','filled'); 
%     if mod(ii,5)==0
%         pp = pp+1; 
%         ax_idx = 4; 
%         plot(ax(ax_idx),y,z)
%         idx(pp) = idx_ii; 
%     end
end
% ax.YTickLabel = num2str(flip(time'));
ax_idx = 2; 
colormap(ax(ax_idx),'jet'); 
cb = colorbar(ax(ax_idx)); 
caxis(ax(ax_idx),[-10,10]); 
set(ax(ax_idx), 'YDir','reverse')
xlabel(ax(ax_idx),'Y pos [mm]');
% ylabel(ax(ax_idx),'Time [#frame]'); 
cb.Label.String = 'd/dt(e_{yy})';
% set(ax(ax_idx), 'XDir','reverse')
%plot(ax(1),st(:,4,:),linspace(Nframe,1,Nframe),'.'); 

% ax_idx = 4; 
% legdata = (squeeze(st(idx,2,end))-squeeze(st(idx_center,2,end)))/resolution; 
% legend(ax(ax_idx),compose("%1.2f",legdata));

% %%
% N = 5; 
% fps = 100; 
% t = 3+3+0.75+3+0.75; 
% signal = repmat(repelem([0 1],N),1,t*fps); 
% 
% time = linspace(0,t,1e3*t); 
% figure; 
% plot(time,signal); 
% xlim([0 0.02]);
% ylim([-0.5 1.5]); 

%% compute repositionning and rescaling of JAI camera 
dir_loc = fullfile(baseDataPath,'rawdata',subject,'calib','2'); 
S = dir(fullfile(dir_loc,sprintf('*cam_0*')));
im_calib = imread(fullfile(S(1).folder,S(1).name));
%
% smallImage = imresize(im_calib, 1/2, 'bilinear');
smallImage = flip(im_calib',2);
figure; 
imshow(smallImage); hold on; 
H = drawcircle; hold on; 
pause; 
c = H.Center; 
r = H.Radius;
drawpoint(gca,'Position',[c(1),c(2)]); hold on; 
%%

%
H = drawline; hold on; 
pause;
%
p1 = H.Position(1,:)'; 
p2 = -H.Position(2,:)'; 
%
dy = -(p2(2)-p1(2));
dx = (p2(1)-p1(1));
alpha = atan2(dx,dy)*180/pi;
%%
dir_loc = fullfile(baseDataPath,'rawdata',subject,'calib','1'); 
S = dir(fullfile(dir_loc,sprintf('*003_undefined_cam_0*')));
im_calib = imread(fullfile(S(1).folder,S(1).name));
figure; 
imshow(im_calib); hold on; 
H = drawline; hold on; 
pause;
p1 = H.Position(:,1)'; 
p2 = H.Position(:,2)'; 
d = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
nbrMM = 24; %%mm from left to right
res = d/nbrMM


%%
BW2 = imfill(masks_saved.cabw_fitr(:,:,end),round([st(idx_center,1,30),st(idx_center,2,30)]));
BW2(round(st(idx_center,1,30)),round(st(idx_center,2,30)))

function flxy = filter_feature_position(flxy)
%FILTER_FEATURE_POSITION 
n = 16; %16; % smoothing parameter (size of gaussian window)
alpha = 10; %10
win = gausswin(n,alpha); win = win/sum(win); % normalized gaussian window

% pad flxy
npad = 4;
flxy = cat(3,repmat(flxy(:,:,1),[1 1 npad]), flxy, repmat(flxy(:,:,end),[1 1 npad]));
flxy = convn(flxy,reshape(win,1,1,n),'same'); % time is the 3rd dim (nframes)
% select center part (remove padding)
flxy = flxy(:,:,npad+1:end-npad);
end


function st = adden(st)
     en = sqrt(st(:,4,:).^2+st(:,5,:).^2+ 2*st(:,6,:).^2);
     st(:,8,:) = en; 
end




