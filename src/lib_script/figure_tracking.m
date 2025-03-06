


baseDataPath = "E:\transfert_SMP"; 
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 

subject = "S03"; 
material = "coating"; 
trial = "001"; 
stereopair = 1; 
phase = "loading";
framejump = 2; 

trial_jj = str2double(trial); 
frameidx = load(fullfile(baseResultPath,'analysis',subject,'frameidx.mat')); frameidx = frameidx.frameidx; 
Fs = 1e3; 
fps = 50; 
idx_frame_start = frameidx{2}.(phase).start;
idx_frame_end = frameidx{2}.(phase).end;

% read images from video file and retreive camera number according to
% stereopair given. 

[cam_first, cam_second, cam_1, cam_2] = import_vid(baseDataPath,...
                                       'subject',subject,...
                                       'material',material,...
                                       'trial',trial,...
                                       'stereopair',stereopair,...
                                       'phase',phase,...
                                       'framejump',framejump,...%);%,...
                                       'idxstart_set',idx_frame_start(trial_jj),... 
                                       'idxend_set',idx_frame_end(trial_jj));
                                   
targetfile = sprintf('ncorr%d.mat',cam_1);
file1 = load(fullfile(baseResultPath,"analysis",...
    subject,...
    material,...
    trial,...
    phase,...
    targetfile)); 
targetfile = sprintf('ncorr%d%d.mat',cam_1,cam_2);
file2 = load(fullfile(baseResultPath,"analysis",...
    subject,...
    material,...
    trial,...
    phase,...
    targetfile));
targetfile = sprintf('ncorr%d.mat',cam_2);
file3 = load(fullfile(baseResultPath,"analysis",...
    subject,...
    material,...
    trial,...
    phase,...
    targetfile)); 
data_ncorr(1) = file1; 
data_ncorr(2) = file2; 
data_ncorr(3) = file3; 
disp('--> reading done')
%%
player(cam_first); 
%%
cam_first = cam_first(:,:,1:end/2+5);  
cam_second = cam_second(:,:,1:end/2+5);  
%% point tracking
% pinit_x = [1328]; 
% pinit_y = [417];
truefps = 50; 
saturationlevel = 500; 
xpos = [600,700,800,900,1000,1100,1200,1300,1400];
ypos = ([385,395,405,415,425,435,445,455,465])+000;
pinit_x = [xpos(4:end),xpos,xpos,xpos,xpos(2:end)]'-400; 
pinit_y = [ypos(4:end)-250,ypos-150,ypos-50,ypos+50,ypos(2:end)+150]'+100;  
pvisu = 2;
pinit = [pinit_x,pinit_y];
sizim = size(cam_first); 
points_match = points_tracking(data_ncorr(2).data_dic_save.displacements(1),pinit,[sizim(1),sizim(2),1]);
points_tracked1 = points_tracking(data_ncorr(1).data_dic_save.displacements,pinit,sizim);
points_tracked2 = points_tracking(data_ncorr(3).data_dic_save.displacements,points_match,sizim);

%format images
%im first 
im = cam_first*4; 
sizim = size(im); 
mask = get_mask(data_ncorr(1).data_dic_save,sizim); 
im_m = im; 
im_m(~mask) = im(~mask).*0.25; 
im_c1 = tracking_feature_im(im_m,points_tracked1,pvisu,...
    'winsize',40,...
    'level',saturationlevel,...
    'pushtop',300+120,...
    'pushright',400); 
im_c1 = im_m; 
cropy = 1:sizim(1);%[1+120:sizim(1)-(300+120)];
cropx = 1:sizim(2);%[1+400:sizim(2)-(200+200)]; 
im_c1 = im_c1(cropy,cropx,:);
points_tracked1 = points_tracked1 + [-cropx(1),-cropy(1)];
%im second
im = cam_second*4;  
sizim = size(im); 
mask = get_mask(data_ncorr(3).data_dic_save,sizim); 
im_m = im; 
im_m(~mask) = im(~mask).*0.25; 
im_c2 = tracking_feature_im(im_m,points_tracked2,pvisu,...
    'winsize',40,...
    'level',saturationlevel,...
    'pushtop',300+50,...
    'pushright',200+150); 
im_c2 = im_m; 
cropy = 1:sizim(1);%1:[1+175:sizim(1)-(300+50)]; 
cropx = 1:sizim(2);%[1+200+150:sizim(2)-(400+200-150)];
im_c2 = im_c2(cropy,cropx,:);
points_tracked2 = points_tracked2 + [size(im_c2,2)-1,0] + [-cropx(1),-cropy(1)];

%anim 
player([im_c1(:,:,1:end),im_c2(:,:,1:end)],...
    'points',[points_tracked1(1,:,:);points_tracked2(1,:,:)],...
    'morepoints',[points_tracked1(1:end,:,:);points_tracked2(1:end,:,:)],...
    'fps',truefps/framejump*0.5,...
    'movie',1);
%%
[cam_first, cam_second, cam_1, cam_2] = import_vid(baseDataPath,...
                                       'subject',subject,...
                                       'material',material,...
                                       'trial',trial,...
                                       'stereopair',1,...
                                       'phase','all',...
                                       'framejump',1);
                                   
[cam_fourth, cam_third, cam_1, cam_2] = import_vid(baseDataPath,...
                                       'subject',subject,...
                                       'material',material,...
                                       'trial',trial,...
                                       'stereopair',2,...
                                       'phase','all',...
                                       'framejump',1);
%%
player(cam_fourth(:,:,1:2:320),'movie',1,'fps',25); 
%%
player(cam_third(:,:,1:2:320),'movie',1,'fps',25); 
%%
player(cam_second(:,:,1:2:320),'movie',1,'fps',25); 
%%
player(cam_first(:,:,1:2:320),'movie',1,'fps',25); 
                                       
%%
function varargout = points_tracking(displ,pinit,sizim,varargin)
[xpos_n_ref, ypos_n_ref] = point_displacement_ref(pinit,displ,sizim);
npoint = size(pinit,1);
N = sizim(3); 
z = zeros(npoint,2,N);
points_ref = z; 
for i = 1:N
    points_ref(:,1,i) = xpos_n_ref(i,:); points_ref(:,2,i) = ypos_n_ref(i,:);
end
if nargout == 1
    varargout{1} = points_ref;
end
end
function [xpos_n, ypos_n] = point_displacement_ref(pos_init,displ,sizim)
N = sizim(3); 
booldispl = 0; 
if N ~= length(displ)
    %it means that it is the second ncorr run, and the first image was
    %compare itself for the first ncorr run only (we take the first image
    %of the first set of image as reference). 
    booldispl = 1; 
end
L = size(pos_init,1);
xpos_n = zeros(N,L); ypos_n = zeros(N,L);
xinit = pos_init(:,1); yinit = pos_init(:,2);
factor = [1216/111,1936/176]; 
for j = 1:L
    xpos = zeros(N,1); ypos = zeros(N,1);
    xpos(1) = xinit(j); ypos(1) = yinit(j);
    if ~booldispl 
        for i = 1:N
            dx = displ(i).plot_u_ref_formatted(round(ypos(1)/factor(1)),round(xpos(1)/factor(2)));
            dy = displ(i).plot_v_ref_formatted(round(ypos(1)/factor(1)),round(xpos(1)/factor(2)));
            xpos(i) = xpos(1)+dx;
            ypos(i) = ypos(1)+dy;
        end
    else 
        for i = 2:N
            dx = displ(i-1).plot_u_ref_formatted(round(ypos(1)/factor(1)),round(xpos(1)/factor(2)));
            dy = displ(i-1).plot_v_ref_formatted(round(ypos(1)/factor(1)),round(xpos(1)/factor(2)));
            xpos(i) = xpos(1)+dx;
            ypos(i) = ypos(1)+dy;
        end
    end
    xpos_n(:,j) = xpos;
    ypos_n(:,j) = ypos;
end
end
function mask = get_mask(data, siz)
    mask = zeros(siz); 
    mask(:,:,1) = imresize(data.displacements(1).roi_cur_formatted.mask,[siz(1),siz(2)]);  
    for ii = 2:siz(3)
        mask(:,:,ii) = imresize(data.displacements(ii-1).roi_cur_formatted.mask,[siz(1),siz(2)]);  
    end
end