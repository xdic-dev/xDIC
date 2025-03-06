function varargout = tracking_feature_im(im, points_tracked, pvisu, varargin)
p = inputParser;
p.addParameter('winsize',20,@isnumeric);
p.addParameter('level',100,@isnumeric);
p.addParameter('pushtop',0,@isnumeric);
p.addParameter('pushright',0,@isnumeric);

p.parse(varargin{:});
winsize = p.Results.winsize;
level = p.Results.level;
pushtop = p.Results.pushtop;
pushright = p.Results.pushright;

border = 10; 
k = 4; 
siz = size(im);
im_c = satur(im,'level',level);
for iframe = 1:siz(3)
    roi = im(round(squeeze(points_tracked(pvisu,2,iframe))-winsize):round(squeeze(points_tracked(pvisu,2,iframe))+winsize), ...
             round(squeeze(points_tracked(pvisu,1,iframe))-winsize):round(squeeze(points_tracked(pvisu,1,iframe))+winsize),iframe); 
    yposwin = siz(1)-(2*k*winsize+(k-1))-border-pushtop:siz(1)-border-pushtop; 

    xposwin = 1+pushright+border:(2*k*winsize+(k-1))+border+pushright+1;

    im_c(yposwin,xposwin,iframe) = imresize(roi,k);
    centerwin = [(yposwin(end)+yposwin(1))/2-1,(xposwin(end)+xposwin(1))/2]; 
    
    Lcross = 75; %pixel
    Lwidth = 7; %pixel
    for ii = round(centerwin(1)-Lcross/2:centerwin(1)+Lcross/2)
        for jj = round(centerwin(2)-Lcross/2:centerwin(2)+Lcross/2)
            if ((ii >= centerwin(1)-Lwidth)&&(ii <= centerwin(1)+Lwidth)||...
                    (jj >= centerwin(2)-Lwidth)&&(jj <= centerwin(2)+Lwidth))
                im_c(ii,jj,iframe) = 255; 
            end
        end
    end
    framewidth = 1; 
    im_c(yposwin(1)-framewidth:yposwin(end)+framewidth,xposwin(1)-framewidth:xposwin(1)+framewidth,:) = 255;
    im_c(yposwin(1)-framewidth:yposwin(end)+framewidth,xposwin(end)-framewidth:xposwin(end)+framewidth,:) = 255;
    im_c(yposwin(1)-framewidth:yposwin(1)+framewidth,xposwin(1)-framewidth:xposwin(end)+framewidth,:) = 255;
    im_c(yposwin(end)-framewidth:yposwin(end)+framewidth,xposwin(1)-framewidth:xposwin(end)+framewidth,:) = 255;
            
end
% % plot of tracking results
% imshow(im_c,[0 im_satur]); hold on; 
% plot(squeeze(points_tracked(:,1,:))',squeeze(points_tracked(:,2,:))','b.-'); hold on;
% plot(squeeze(points_tracked(:,1,iframe)),squeeze(points_tracked(:,2,iframe)),'yo','MarkerSize',5);
if nargout == 1
    varargout{1} = im_c;
end
end