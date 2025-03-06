function varargout = player(imarray,varargin)
%PLAYER Image sequence visualization.
%Tool to observe an image sequence with several elements superimposed,
%which are passed as arguments (see inputs).
%
% Syntax: h = player(imarray,varargin)
%
% Example: h = player(imarray,'contour',caxy,'points',flxy);
%
% Inputs:
%   imarray     array of pictures                              dim = (m,n,npic)
%   contour     [opt] xy coord of the contour of contact area  dim = (npts,2,npic)
%                          OR bw image of the contact area     dim = size(imarray(:,:,1))
%   ellipse     [opt] parameters of ellipses fitted on CA      dim = (1,5,npic)
%   points      [opt] xy coord of optical flow points          dim = (npts,2,npic)
%   morepoints  [opt] xy coord of other points                 dim = (npts,2,npic)
%   vectors     [opt] xy coord and uv components of vectors    dim = (npts,4,npic)
%   strain      [opt] array to display as in strain_show.m     dim = (see below)
%   flipcmap    [opt] flip color map (as in strain_show.m)     default: false
%   subreg      [opt] subregion to search for CA [x y w h]     default: [0 0 0 0]
%   npts        [opt] number of points in CA boundary          default: 50
%   title       [opt] title of the player window               default: 'SMP Video Player'
%   labels      [opt] display stimulation direction labels     default: false
%   fps         [opt] playback speed (frames per second)       default: 100 fps
%   movie       [opt] if true, enable movie recording          default: false
%   moviesize   [opt] x and y dimensions of the movie          default: [800,640]
%                        note: x/y ratio must be 1.25 (as the images)
%                        examples: [800,640] (default) ; [1280,1024]
%   allon       [opt] if true, display opt items at startup    default: false
%   winsize     [opt] dims of the player window [x y w h]      default: []
%
% Outputs:
%   h           [opt] figure handle (useful to delete GUI, e.g.)
%
% Strain values displayed in the player
%   The 'strain' optional input allows the user to pass an array of values
%   which will be displayed as in strain_show.m, i.e. using colors. This
%   array must contain the x and y position of the triangle centers along
%   with the values to display in this triangles. Therefore, the array has
%   3 components: x, y and the data to display.
%   The player accepts the follwing input dimensions:
%       (ntris,3,npic)    -->   will be displayed as is
%       (ntris,3,<npic)   -->   the last frame will be repeated until the
%                               3rd dim of the strain array matches npic.
%
% Player functionalities and keyboard shortcuts:
%   [Space]                 Play/Pause
%   [Ctrl + Space]          Stop
%   [l]                     Toggle loop mode
%
%   [Right]                 Go to next frame (= jump 1 frame forward)
%   [Ctrl + Right]          Jump 10 frames forward
%   [Ctrl + Shift + Right]  Jump 100 frames forward
%
%   [Left]                  Go to previous frame (= jump 1 frame backward)
%   [Ctrl + Left]           Jump 10 frames backward
%   [Ctrl + Shift + Left]   Jump 100 frames backward
%
%   [Up arrow]              FPS x2 (max is 200 fps)
%   [Down arrow]            FPS x0.5
%
%   [k]                     Increase Skip of 1 frame
%   [Shift + k]             Decrease Skip of 1 frame (min is 0)
%
%   [v]                     Toggle movie mode and reopen file
%   [s]                     Snapshot of current frame (uncompressed TIFF)
%   [r]                     Reset settings (Current frame, FPS and Skip)
%   [t]                     Toggle stimulation direction text labels
%
%   [c]                     Toggle contour display
%   [e]                     Toggle ellipse display
%   [p]                     Toggle points display
%   [m]                     Toggle morepoints display
%   [a]                     Toggle strain-like display
%   [b]                     Toggle subregion display
%
%   [Escape or q]           Quit player
%
% Editable fields:
%   Current frame: can be edited to go to any frame directly.
%
%   FPS: indicates playback speed. The default is 100 frames per second.
%
%   Skip: indicates how many frames must be skipped between two displayed
%     frames. The default is 0, i.e. display every frame.
%
% Movie recording mode:
%   When the 'movie' parameter is true, the player is in movie recording
%   mode. The movie will capture every elements displayed on the figure,
%   such as superimposed points. Note also that the movie quality will
%   depend on the actual size of the figure.
%
%   When movie mode is enabled, the movie starts recording when the user
%   presses the Play button. Once the last frame is reached, a pop-up
%   window will show up, asking for the output file path. The movie is not
%   recorded if the user hits the Cancel button on the pop-up window. The
%   movie is encoded as an MP4 video using the H.264 encoder.
%
%   The movie encoding is done through imarray2mp4.m, using avconv.
%   To use this function, you must have the avconv executable in your
%   system PATH environment variable.
%
% Note on stimulation directions:
%   The images displayed in this player are flipped vertically, i.e. along
%   the y-axis. This is due to the way the Matlab 'imshow' function works.
%   To clarify stimulation directions, text labels have been added in the
%   player. The x-direction is DP (distal-proximal) and the y-direction is
%   RU (radial-ulnar). Radial is towards the thumb side of the hand and
%   Proximal is towards the wrist.

% Parse inputs
p = inputParser;
p.addParameter('contour',[]);
p.addParameter('ellipse',[]);
p.addParameter('points',[]);
p.addParameter('morepoints',[]);
p.addParameter('vectors',[]);
p.addParameter('strain',[]);
p.addParameter('flipcmap',true);
p.addParameter('subreg',[]);
p.addParameter('npts',50);
p.addParameter('title','SMP Video Player');
p.addParameter('labels',false);
p.addParameter('fps',100);
p.addParameter('movie',false);
p.addParameter('moviesize',[2*2560/3,2048/3]); % Felicien: 640 doit etre change en 590 selon le PC!!!
p.addParameter('allon',false);
p.addParameter('winsize',[]);
p.addParameter('imposedStLimit',[]);
p.addParameter('imposedImLimit',[]);

p.parse(varargin{:});
contour    = p.Results.contour;
ellipse    = p.Results.ellipse;
points     = p.Results.points;
morepoints = p.Results.morepoints;
vectors    = p.Results.vectors;
strain     = p.Results.strain;
flipcmap   = p.Results.flipcmap;
subreg     = p.Results.subreg;
npts       = p.Results.npts;
title      = p.Results.title;
labels     = p.Results.labels;
fps        = p.Results.fps;
movieon    = p.Results.movie;
moviesize  = p.Results.moviesize;
allon      = p.Results.allon;
winsize    = p.Results.winsize;
imposedStLimit = p.Results.imposedStLimit;
imposedImLimit = p.Results.imposedImLimit;

if(~isempty(contour)),  contour_visible = 'on';
else  contour_visible = 'off'; end
if(~isempty(ellipse)),  ellipse_visible = 'on';
else  ellipse_visible = 'off'; end
if(~isempty(points)),  points_visible = 'on';
else  points_visible = 'off'; end
if(~isempty(morepoints)),  morepoints_visible = 'on';
else  morepoints_visible = 'off'; end
if(~isempty(vectors)),  vectors_visible = 'on';
else  vectors_visible = 'off'; end
if(~isempty(strain)), strain_visible = 'on';
else  strain_visible = 'off'; end
if(~isempty(subreg)),  subreg_visible = 'on';
else  subreg_visible = 'off'; end
if(isempty(subreg)),  subreg = [0,0,eps,eps]; end

% Output array
mOutputArgs = {};

% Movie length
npic = size(imarray,3);

% If contour is a bw image, extract xy coord of CA boundary
if(isinteger(contour) || islogical(contour))
  sarr = size(imarray);
  scont = size(contour);
  % Also works if contour contains a single image
  if(length(scont) == 2)
    scont(3) = 1;
  end
  xyconv = zeros(npts,2,scont(3));
  for jj = 1:scont(3)
    if(isa(contour,'uint32')) % if bw is packed
      loc = bwunpack(contour(:,:,jj));
    else
      loc = contour(:,:,jj);
    end
    loc = edge(imresize(loc,sarr(1:2)),'zerocross');
    clear xys;
    [xys(1), xys(2)] = find(loc,1);
    % xyconv(:,[2 1],jj) = tab2convex(xys,npts);
    xy = bwtraceboundary(loc,[xys(1), xys(2)],'N');
    xyconv(:,[2 1],jj) = interp1(linspace(0,1,length(xy)),xy,...
      linspace(0,1,npts));
  end
  contour = xyconv;
end

% Process strain input, if needed
if(~isempty(strain))
  % Repeat last frame to match imarray length
  siz_st = size(strain,3);
  if(siz_st < npic)
    strain = cat(3,strain,repmat(strain(:,:,end),1,1,npic-siz_st));
  end
  
  % Extract x, y and data and compute triangulation
  stX = squeeze(double(strain(:,1,:)));
  stY = squeeze(double(strain(:,2,:)));
  stDataRaw = squeeze(double(strain(:,3,:)));
  stTri = delaunay(stX(:,1),stY(:,1));
  
  % Define number of color levels and color map
  stColorLvls = 256;
  if(flipcmap)
    stCmap = [1 1 1; flipud(jet(stColorLvls))]; % NaN -> 0, displayed in white
  else
    stCmap = [1 1 1; jet(stColorLvls)]; % NaN -> 0, displayed in white
  end
  
  % Normalize strain data to max on 'stColorLvls' color level
  if(isempty(imposedStLimit))
    stLimit = max(abs(prctile(stDataRaw(:),[1,99])));
  else
    stLimit=imposedStLimit;
  end
  stData = ceil((stDataRaw+stLimit)/(2*stLimit)*stColorLvls);
  stData(stData > stColorLvls) = stColorLvls;
  stData(stData < 1) = 1;
  stData(isnan(stData)) = 0; % NaN -> 0, displayed in white
  stData = stData + 1;
end

% Parameters
fr = 1;
pauseRequest = true;
loopStatus = false;
initfps = fps;
skip = 0;

% UI objects
ui.fh = figure('Name',title,'NumberTitle','off', ...
  'Color',get(0,'defaultuicontrolbackgroundcolor'));
if(~isempty(winsize)), set(ui.fh,'Position',winsize); end
set(ui.fh,'Toolbar','Figure'); % prevent toolbar removing by uicontrol
set(ui.fh,'KeyPressFcn',@keypresscb); % register key press listener
ui.ax = gca;

% Adjust figure size for movie (assuming FullHD screen (1920x1080))
% Note: default OuterPosition prop of axis = [0.1300 0.1100 0.7750 0.8150]
if(movieon)
  screenw = 1920;
  screenh = 1080;
  
  figw = round(moviesize(1)/0.7750);
  figh = round(moviesize(2)/0.8150);
  figl = max(floor((screenw-figw)/2),0);
  figb = max(floor((screenh-figh)/2),0);
  
  axw = moviesize(1);
  axh = moviesize(2);
  axl = round(0.1300*figw);
  axb = round(0.1100*figh);
  
  set(ui.fh,'Units','pixels','Position',[figl,figb,figw,figh]);
  set(ui.ax,'Units','pixels','Position',[axl,axb,axw,axh]);
end

ui.prev100Button = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','<<<','callback', ...
  {@changeimcb,-100},'Position',[.13 .06 .09 .04]);
ui.prev10Button = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','<<','callback', ...
  {@changeimcb,-10},'Position',[.22 .06 .09 .04]);
ui.prevButton = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','<','callback', ...
  {@changeimcb,-1},'Position',[.31 .06 .09 .04]);

ui.playButton = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','Play','callback', ...
  @playcb,'Position',[.4 .06 .078 .04]);
ui.stopButton = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','Stop','callback', ...
  @stopcb,'Position',[.478 .06 .079 .04]);
ui.loopButton = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','Loop OFF','callback', ...
  @loopcb,'Position',[.557 .06 .078 .04]);

ui.nextButton = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','>','callback', ...
  {@changeimcb,1},'Position',[.635 .06 .09 .04]);
ui.next10Button = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','>>','callback', ...
  {@changeimcb,10},'Position',[.725 .06 .09 .04]);
ui.next100Button = uicontrol(ui.fh,'Style','pushbutton', ...
  'Units','normalized','String','>>>','callback', ...
  {@changeimcb,100},'Position',[.815 .06 .09 .04]);

ui.currimStr = uicontrol(ui.fh,'Style','edit', ...
  'String','','callback',@changeimcb, ...
  'Units','normalized','Position',[.3 .94 .1 .04]);
ui.currimLabel = uicontrol(ui.fh,'Style','text','String','Current frame:', ...
  'Units','normalized','Position',[.15 .935 .15 .04]);
ui.fpsStr = uicontrol(ui.fh,'Style','edit', ...
  'String',num2str(fps),'callback',@fpscb, ...
  'Units','normalized','Position',[.535 .94 .1 .04]);
ui.fpsLabel = uicontrol(ui.fh,'Style','text','String','FPS:', ...
  'Units','normalized','Position',[.425 .935 .1 .04]);
ui.skipStr = uicontrol(ui.fh,'Style','edit', ...
  'String',num2str(skip),'callback',@skipcb, ...
  'Units','normalized','Position',[.75 .94 .1 .04]);
ui.skipLabel = uicontrol(ui.fh,'Style','text','String','Skip:', ...
  'Units','normalized','Position',[.65 .935 .1 .04]);

ui.contourCheck = uicontrol(ui.fh,'Style','check','value',0, ...
  'String','contour','callback',@contourCheckcb, ...
  'Units','normalized','Position',[.14 .01 .2 .04], ...
  'Enable',contour_visible);
ui.ellipseCheck = uicontrol(ui.fh,'Style','check','value',0, ...
  'String','ellipse','callback',@ellipseCheckcb, ...
  'Units','normalized','Position',[.27 .01 .2 .04], ...
  'Enable',ellipse_visible);
ui.pointsCheck = uicontrol(ui.fh,'Style','check','value',0, ...
  'String','points','callback',@pointsCheckcb, ...
  'Units','normalized','Position',[.40 .01 .2 .04], ...
  'Enable',points_visible);
ui.morepointsCheck = uicontrol(ui.fh,'Style','check','value',0, ...
  'String','more points','callback',@morepointsCheckcb, ...
  'Units','normalized','Position',[.52 .01 .2 .04], ...
  'Enable',morepoints_visible);
ui.vectorsCheck = uicontrol(ui.fh,'Style','check','value',0, ...
  'String','vectors','callback',@vectorsCheckcb, ...
  'Units','normalized','Position',[.67 .01 .2 .04], ...
  'Enable',vectors_visible);
ui.strainCheck = uicontrol(ui.fh,'Style','check','value',0, ...
  'String','strain','callback',@strainCheckcb, ...
  'Units','normalized','Position',[.79 .01 .2 .04], ...
  'Enable',strain_visible);
% ui.subregCheck = uicontrol(ui.fh,'Style','check','value',0, ...
%   'String','subregion','callback',@subregCheckcb, ...
%   'Units','normalized','Position',[.79 .01 .2 .04], ...
%   'Enable',subreg_visible);

% Plot graphic objects
%display range is set to min and max values of the image
himg  = imshow(imarray(:,:,end),imposedImLimit,'Parent',ui.ax);
set(ui.ax,'nextplot','add');

% if ~isempty(points)
%     plot(squeeze(points(:,1,:))',squeeze(points(:,2,:))','Color',"#FFFF00",'LineWidth',3); %,"#FFFF00"
% end
% if ~isempty(morepoints)
%     plot(squeeze(morepoints(:,1,:))',squeeze(morepoints(:,2,:))','Color',"#FFFF00",'LineWidth',3); %,"#FFFF00"
% end
hpts  = plot(nan,nan,'r+','Linewidth',1,'Parent',ui.ax,...
    'markersize',1,'LineWidth',0.5,'MarkerFaceColor','k');%,'LineWidth',5);%15
hmpts = plot(nan,nan,'b+','Linewidth',1,'Parent',ui.ax,...
    'markersize',1,'LineWidth',0.5,'MarkerFaceColor','k');
hvec = quiver(nan,nan,nan,nan,0,'c','Linewidth',1,'Parent',ui.ax);
hcont = plot(nan,nan,'m','Linewidth',2,'Parent',ui.ax);
hell  = plot(nan,nan,'y','Linewidth',2,'Parent',ui.ax);
hrect = rectangle('Position',subreg,'EdgeColor',[0,0.75,0],...
  'Linewidth',3,'Parent',ui.ax);
set([hpts,hmpts,hvec,hcont,hell,hrect],'Visible','off');

if(~isempty(strain))
  % Create patch object
  hst = patch('Faces',stTri,'Vertices',nan(size(stTri,1),2),...
    'FaceColor','interp','FaceVertexCData',nan(size(stTri,1),1),...
    'FaceVertexAlphaData',ones(size(stTri,1),1),...
    'AlphaDataMapping','none',...
    'EdgeColor','none','FaceAlpha','interp','Parent',ui.ax);
  set(hst,'Visible','off');
  
  % Create custom axes for custom colorbar
  ui.axStCmap = axes('Parent',ui.fh,'Visible','off');
  uistack(ui.axStCmap,'top');
  if(flipcmap)
    colormap(ui.axStCmap,flipud(jet(stColorLvls)));
  else
    colormap(ui.axStCmap,jet(stColorLvls));
  end
  caxis(ui.axStCmap,[-stLimit,stLimit]);
end

% Display stimulation direction text labels
sizX = size(imarray,2);
sizY = size(imarray,1);
htxtU = text(sizX/2,25,'Ulnar','FontSize',12,'Color',[1,1,1],...
  'HorizontalAlignment','center','Parent',ui.ax);
htxtD = text(25,sizY/2,'Distal','FontSize',12,'Color',[1,1,1],...
  'HorizontalAlignment','center','Rotation',90,'Parent',ui.ax);
htxtR = text(sizX/2,sizY-34,'Radial','FontSize',12,'Color',[1,1,1],...
  'HorizontalAlignment','center','Parent',ui.ax);
htxtP = text(sizX-40,sizY/2,'Proximal','FontSize',12,'Color',[1,1,1],...
  'HorizontalAlignment','center','Rotation',90,'Parent',ui.ax);
if(~labels), set([htxtU,htxtD,htxtR,htxtP],'Visible','off'); end

% Start player
changeimcb(0,0,0);
if(allon)
  toggleContour();
  toggleEllipse();
  togglePoints();
  toggleMorepoints();
  toggleVectors();
  toggleStrain();
  toggleSubregion();
end

% Return output if necessary
mOutputArgs{1} = ui.fh;
if(nargout > 0)
  [varargout{1:nargout}] = mOutputArgs{:};
end

%------------------------------------------------------------------------
  function playcb(varargin)
    if(pauseRequest)
      pauseRequest = false;
      set(ui.playButton,'String','Pause');
      if(fr >= npic)
        fr = 1;
        set(ui.currimStr,'String',num2str(fr));
      end
      if(movieon)
        F(length(fr:npic)) = struct('cdata',[],'colormap',[]);
      end
      for ii = 1:length(fr:npic)
        if(movieon), F(ii) = getframe(ui.ax); end
        changeimcb(0,0,1+skip); pause(1/fps);
        if(pauseRequest), break; end
      end
      if(movieon && ~pauseRequest)
        savemovie(F);
      end
      if(fr >= npic)
        pauseRequest = true;
        if(loopStatus)
          playcb(); % recursive call, must be the last instruction!
        else
          set(ui.playButton,'String','Play');
        end
      end
    else
      pauseRequest = true;
      set(ui.playButton,'String','Play');
    end
  end

%------------------------------------------------------------------------
  function stopcb(varargin)
    pauseRequest = true;
    set(ui.playButton,'String','Play');
    fr = 1;
    set(ui.currimStr,'String',num2str(fr));
    changeimcb(0,0,0);
  end

%------------------------------------------------------------------------
  function loopcb(varargin)
    if(loopStatus)
      loopStatus = false;
      set(ui.loopButton,'String','Loop OFF');
    else
      loopStatus = true;
      set(ui.loopButton,'String','Loop ON');
    end
  end

%------------------------------------------------------------------------
  function fpscb(varargin)
    fps = str2double(get(ui.fpsStr,'String'));
    fps = min(fps,200);
    set(ui.fpsStr,'String',num2str(fps));
  end

%------------------------------------------------------------------------
  function modiffps(inc)
    fps = fps * inc;
    fps = min(fps,200);
    set(ui.fpsStr,'String',num2str(fps));
  end

%------------------------------------------------------------------------
  function skipcb(varargin)
    skip = str2double(get(ui.skipStr,'String'));
    skip = floor(skip);
    skip = max(skip,0);
    set(ui.skipStr,'String',num2str(skip));
  end

  %------------------------------------------------------------------------
  function modifskip(inc)
    skip = floor(skip + inc);
    skip = max(skip,0);
    set(ui.skipStr,'String',num2str(skip));
  end

  %------------------------------------------------------------------------
  function contourCheckcb(varargin)
    changeimcb();
    if(get(ui.contourCheck,'value'))
      set(hcont,'Visible','on');
    else
      set(hcont,'Visible','off');
    end
  end

%------------------------------------------------------------------------
  function toggleContour()
    if(strcmp(contour_visible,'on'))
      set(ui.contourCheck,'value',~get(ui.contourCheck,'value'));
      contourCheckcb();
    end
  end

%------------------------------------------------------------------------
  function ellipseCheckcb(varargin)
    changeimcb();
    if(get(ui.ellipseCheck,'value'))
      set(hell,'Visible','on');
    else
      set(hell,'Visible','off');
    end
  end

%------------------------------------------------------------------------
  function toggleEllipse()
    if(strcmp(ellipse_visible,'on'))
      set(ui.ellipseCheck,'value',~get(ui.ellipseCheck,'value'));
      ellipseCheckcb();
    end
  end

%------------------------------------------------------------------------
  function pointsCheckcb(varargin)
    changeimcb();
    if(get(ui.pointsCheck,'value'))
      set(hpts,'Visible','on');
    else
      set(hpts,'Visible','off');
    end
  end

%------------------------------------------------------------------------
  function togglePoints()
    if(strcmp(points_visible,'on'))
      set(ui.pointsCheck,'value',~get(ui.pointsCheck,'value'));
      pointsCheckcb();
    end
  end

%------------------------------------------------------------------------
  function morepointsCheckcb(varargin)
    changeimcb();
    if(get(ui.morepointsCheck,'value'))
      set(hmpts,'Visible','on');
    else
      set(hmpts,'Visible','off');
    end
  end

%------------------------------------------------------------------------
  function toggleMorepoints()
    if(strcmp(morepoints_visible,'on'))
      set(ui.morepointsCheck,'value',~get(ui.morepointsCheck,'value'));
      morepointsCheckcb();
    end
  end
%------------------------------------------------------------------------
  function vectorsCheckcb(varargin)
    changeimcb();
    if(get(ui.vectorsCheck,'value'))
      set(hvec,'Visible','on');
    else
      set(hvec,'Visible','off');
    end
  end

%------------------------------------------------------------------------
  function toggleVectors()
    if(strcmp(vectors_visible,'on'))
      set(ui.vectorsCheck,'value',~get(ui.vectorsCheck,'value'));
      vectorsCheckcb();
    end
  end

%------------------------------------------------------------------------
  function strainCheckcb(varargin)
    changeimcb();
    if(get(ui.strainCheck,'value'))
      set(hst,'Visible','on');
      if(isempty(imposedStLimit))
        colorbar(ui.axStCmap,...
          'AxisLocation','out',...
          'Location','manual',...
          'Position',[0.07,0.1100,0.04,0.8150],...
          'Units','normalized',...
          'Visible','on');
      else
        colorbar(ui.axStCmap,...
          'AxisLocation','out',...
          'Location','manual',...
          'Position',[0.07,0.1100,0.04,0.8150],...
          'Units','normalized',...
          'Limits',[-imposedStLimit imposedStLimit],...
          'Visible','on');
      end
    else
      set(hst,'Visible','off');
      colorbar(ui.axStCmap,'off');
    end
  end

%------------------------------------------------------------------------
  function toggleStrain()
    if(strcmp(strain_visible,'on'))
      set(ui.strainCheck,'value',~get(ui.strainCheck,'value'));
      strainCheckcb();
    end
  end

%------------------------------------------------------------------------
  function subregCheckcb(varargin)
    if(get(ui.subregCheck,'value'))
      set(hrect,'Visible','on');
    else
      set(hrect,'Visible','off');
    end
  end

%------------------------------------------------------------------------
  function toggleSubregion()
    if(strcmp(subreg_visible,'on'))
      set(ui.subregCheck,'value',~get(ui.subregCheck,'value'));
      subregCheckcb();
    end
  end

%------------------------------------------------------------------------
  function toggleLabels()
    if(labels)
      labels = false;
      set([htxtU,htxtD,htxtR,htxtP],'Visible','off');
    else
      labels = true;
      set([htxtU,htxtD,htxtR,htxtP],'Visible','on');
    end
  end

%------------------------------------------------------------------------
  function changeimcb(varargin)
    % Compute current frame number
    if(length(varargin) < 3)
      fr = floor(str2double(get(ui.currimStr,'String')));
    else
      inc = varargin{3};
      fr = fr + inc;
    end
    if(fr < 1), fr = 1; end
    if(fr > npic), fr = npic; end
    set(ui.currimStr,'String',num2str(fr));
    
    % Update graphic objects
    set(himg,'cdata',imarray(:,:,fr));
    if(get(ui.contourCheck,'value'))
      set(hcont,'xdata',contour(:,1,fr),'ydata',contour(:,2,fr));
    end
    if(get(ui.ellipseCheck,'value'))
      params = ellipse(1,:,fr);
      t_ell = linspace(0,pi*2);
      x_ell = params(3)*cos(t_ell);
      y_ell = params(4)*sin(t_ell);
      nx = x_ell*cos(params(5)) - y_ell*sin(params(5)) + params(1);
      ny = x_ell*sin(params(5)) + y_ell*cos(params(5)) + params(2);
      set(hell,'xdata',nx,'ydata',ny);
    end
    if(get(ui.pointsCheck,'value'))
      if(fr <= size(points,3))
        set(hpts,'xdata',points(:,1,fr),'ydata',points(:,2,fr));
      end
    end
    if(get(ui.morepointsCheck,'value'))
      if(fr <= size(morepoints,3))
        set(hmpts,'xdata',morepoints(:,1,fr),'ydata',morepoints(:,2,fr));
      end
    end
    if(get(ui.vectorsCheck,'value'))
      if(fr <= size(vectors,3))
        set(hvec,'xdata',vectors(:,1,fr),'ydata',vectors(:,2,fr),'udata',vectors(:,3,fr),'vdata',vectors(:,4,fr));
      end
    end
    if(get(ui.strainCheck,'value'))
      set(hst,'Vertices',[stX(:,fr),stY(:,fr)],...
        'FaceVertexCData',stCmap(stData(:,fr),:),...
        'FaceVertexAlphaData',(stData(:,fr)>1)*.4);
    end
  end

%------------------------------------------------------------------------
  function saveimage(varargin)
    pauseRequest = true;
    set(ui.playButton,'String','Play');
    [filename,pathname] = uiputfile({'*.tiff','Images (*.tiff)'; ...
      '*.*','All Files (*.*)'},'Save image as');
    if(~isequal(filename,0) && ~isequal(pathname,0))
      imwrite(imscale(imarray(:,:,fr)),fullfile(pathname,filename), ...
        'tiff','Compression','none');
    end
  end

%------------------------------------------------------------------------
  function savemovie(F)
    [filename,pathname] = uiputfile({'*.mp4','Movies (*.mp4)'; ...
      '*.*','All Files (*.*)'},'Save movie as');
    if(~isequal(filename,0) && ~isequal(pathname,0))
%       npicmov = length(F);
%       moviearray = uint8(zeros(moviesize(2),moviesize(1),3,npicmov));
%       for ii = 1:npicmov
%         tmp = F(ii).cdata; % dims of captured frames might not be exact
%         moviearray(:,:,:,ii) = tmp(1:moviesize(2),1:moviesize(1),:);
%       end
      % imarray2mp4(moviearray,fullfile(pathname,filename),false,false);
      v = VideoWriter(fullfile(pathname,filename),'MPEG-4');
      v.Quality = 100;
      v.FrameRate = fps;
%       v.CompressionRatio = 2; 
      open(v);
      
      for ii=1:length(F)
        writeVideo(v, F(ii));
      end
      
      close(v)
    end
  end

%------------------------------------------------------------------------
  function resetplayer(varargin)
    fr = 1;
    pauseRequest = true;
    loopStatus = false;
    fps = initfps;
    skip = 0;
    set(ui.currimStr,'String',num2str(fr));
    set(ui.fpsStr,'String',num2str(fps));
    set(ui.skipStr,'String',num2str(skip));
    set(ui.playButton,'String','Play');
    set(ui.loopButton,'String','Loop OFF');
    changeimcb(0,0,0);
  end

%------------------------------------------------------------------------
% function reopenplayer(varargin)
%   ison = @(x) strcmp(x,'on');
%   args_bool = [ison(contour_visible),ison(ellipse_visible),...
%     ison(points_visible),ison(morepoints_visible),...
%     ison(subreg_visible)]';
%   args = {'contour',contour;'ellipse',ellipse;'points',points;...
%     'morepoints',morepoints;'subreg',subreg};
%   args = args(args_bool,:);
%   args = reshape(args',numel(args),1);
%   player(imarray,args{:},'npts',npts,'title',title,...
%     'fps',initfps,'movie',movieon);
%   close(ui.fh);
% end

%------------------------------------------------------------------------
% Same function as the one commented above, but simpler to understand
  function reopenplayer(varargin)
    if(strcmp(subreg_visible,'on'))
      player(imarray,'contour',contour,'ellipse',ellipse,...
        'points',points,'morepoints',morepoints,...
        'strain',strain,'flipcmap',flipcmap,'subreg',subreg,...
        'npts',npts,'title',title,'labels',labels,'fps',initfps,...
        'movie',movieon,'allon',allon,'winsize',winsize);
    else
      player(imarray,'contour',contour,'ellipse',ellipse,...
        'points',points,'morepoints',morepoints,...
        'strain',strain,'flipcmap',flipcmap,...
        'npts',npts,'title',title,'labels',labels,'fps',initfps,...
        'movie',movieon,'allon',allon,'winsize',winsize);
    end
    close(ui.fh);
  end

%------------------------------------------------------------------------
% function file = pickfile(varargin)
%   [fname,dirpath] = uigetfile('*.bin','Select a camera binary file');
%   fpath = fullfile(dirpath,fname);
%   if(exist(fpath,'file'))
%     file = fpath;
%   else
%     file = [];
%   end
% end

%------------------------------------------------------------------------
  function keypresscb(~,eventdata)
    key = eventdata.Key; % character corresponding to pressed key
    mod = eventdata.Modifier; % cell array, e.g. {'shift','control'}
    shift_bool = any(strcmp('shift',mod));
    ctrl_bool = any(strcmp('control',mod));
    switch(key)
      case 'space'
        if(~ctrl_bool)     % no mod
          playcb();
        elseif(ctrl_bool)  % control
          stopcb();
        end
      case 'l'
        loopcb();
      case 'rightarrow'
        if(~ctrl_bool && ~shift_bool)     % no mod
          changeimcb(0,0,1);
        elseif(ctrl_bool && ~shift_bool)  % control
          changeimcb(0,0,10);
        elseif(ctrl_bool && shift_bool)   % control + shift
          changeimcb(0,0,100);
        end
      case 'leftarrow'
        if(~ctrl_bool && ~shift_bool)     % no mod
          changeimcb(0,0,-1);
        elseif(ctrl_bool && ~shift_bool)  % control
          changeimcb(0,0,-10);
        elseif(ctrl_bool && shift_bool)   % control + shift
          changeimcb(0,0,-100);
        end
      case 'uparrow'
        modiffps(2);
      case 'downarrow'
        modiffps(0.5);
      case 'k'
        if(~shift_bool)     % no mod
          modifskip(1);
        elseif(shift_bool)  % shift
          modifskip(-1);
        end
        % Disabled because not useful.
        % case 'o'
        %   file = pickfile();
        %   if(~isempty(file))
        %     imarray = readbin(file);
        %     reopenplayer();
        %   end
      case 'v'
        if(movieon)
          movieon = false;
        else
          movieon = true;
        end
        reopenplayer();
      case 's'
        saveimage();
      case 'r'
        resetplayer();
      case 't'
        toggleLabels();
      case 'c'
        toggleContour();
      case 'e'
        toggleEllipse();
      case 'p'
        togglePoints();
      case 'm'
        toggleMorepoints();
      case 'a'
        toggleStrain();
      case 'b'
        toggleSubregion();
      case {'q','escape'}
        close(ui.fh);
    end
  end

end
