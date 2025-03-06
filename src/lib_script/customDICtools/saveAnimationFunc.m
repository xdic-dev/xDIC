function saveAnimationFunc(varargin)
%% Parse input arguments
p = inputParser;
p.addParameter('basePath',[]);
p.addParameter('filename',[]);
p.addParameter('format','eps');
p.addParameter('res',120);
p.addParameter('figopt',[]);
p.addParameter('fig',[]);
p.addParameter('movie',0);
p.addParameter('fps',50);
p.addParameter('IDimage',[]);
p.addParameter('target_obj',[]);


p.parse(varargin{:});
basePath = p.Results.basePath;
filename = p.Results.filename;
format = p.Results.format;
res = p.Results.res;
figopt = p.Results.figopt;
figon= p.Results.fig;
movieon = p.Results.movie;
fps = p.Results.fps;
IDimage = p.Results.IDimage;
target_obj = p.Results.target_obj;

%% Setup figure and axes for the recording
hf = gcf; 
ax = gca; 
% set figure visible for fullscreen plot
set(hf,'visible',1);
% set figure fullscreen
set(hf,'WindowState','maximized'); 
% % set the same aspect ratio for x and y axis 
% ax.DataAspectRatio = [1,1,1]; 
% %figure size
% set(ax,'Units', 'centimeters');
% width = 40; 
% OutPos = get(ax,'OuterPosition'); ratio = OutPos(4)/OutPos(3); 
% set(ax,'OuterPosition',[5 2 width width*ratio]);
% set(hf, 'Position', get(0, 'Screensize'));
% % set no title 
% title_ax = get(ax,'Title');
% title_ax = title_ax.String;
% ax.Title.String = []; 
% % set no axes 
%set(ax,'visible','off');
% % set no colobar
colorbar(ax,'hide'); 

% set(hf,'Color','white'); 

% ax_limits(1:2) = get(gca,'XLim');
% ax_limits(3:4) = get(gca,'YLim');
% fig_limits = axescoord2figurecoord(get(gca,'XLim'),get(gca,'YLim'));

%% Slider to change current image at the screen to save
jSlider=hf.UserData.anim8.sliderHandles{1};

% Filename 
%Make sure folder exist
if ~exist(basePath,'dir')
    mkdir(basePath); 
end

%IDimage
numberOfFrames = numel(hf.UserData.anim8.animStruct.Time);
if isempty(IDimage)
    IDimage = 1:numberOfFrames; 
end

imname = sprintf('%s',filename);
editformat=['-',format];
editimRes=sprintf('-r%d',res);
Q = {basePath,imname,editformat,editimRes};
if ~isempty(figopt)
    Q{end+1} = figopt; 
end
completeFileName=fullfile(Q{1},Q{2});

%% 
% if ~exist(fullfile(basePath,'fig'),'dir')
%     mkdir(fullfile(basePath,'fig')); 
%     mkdir(fullfile(basePath,'vid')); 
% end
%% Recording
if movieon == 1
    %Write video 
    videofilename = sprintf('%s.mp4',filename);
    writerObj = VideoWriter(fullfile(basePath,videofilename),'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = fps;

    % open the video writer
    open(writerObj);
end
idx = 1; 
for ii=1:numberOfFrames
    set(jSlider,'Value',ii);
    if strcmp(target_obj,'fig')
        target_obj = gcf; 
    elseif strcmp(target_obj,'ax')
        target_obj = gca;
    end
    if length(IDimage)>=idx
        if figon == 1 && IDimage(idx) == ii 
            fileNameNow=fullfile(basePath,sprintf('%s_%d',filename,ii));
    %         Q{1}=fileNameNow;
            %myexport_fig(Q{:});
%             export_fig(fileNameNow, '-png', '-eps', '-transparent');
            exportgraphics(target_obj,sprintf('%s.%s',fileNameNow,format),...
                'BackgroundColor','none',...
                'ContentType','vector'); 
            idx = idx+1; 
        end
    end
    if movieon == 1
        % write the frames to the video
        thisFrame = getframe(target_obj); 
        % Write this frame out to a new video file.
        writeVideo(writerObj, thisFrame);
    end
end

if movieon == 1
    % close the writer object
    close(writerObj);
    fprintf('movie saved\n'); 
end
% target_obj.Title.String = sprintf('%s\n%s\n%s',title_ax{1},title_ax{2},title_ax{3}); 
% set(target_obj,'visible','on');
% set(hf,'visible',1);
% colorbar(target_obj,'on'); 
end