 
baseDataPath = "E:\transfert_SMP"; 
datafilepath = fullfile(baseDataPath,"rawdata","S03","speckles","coating","robot");
data = import_robotdata_subject(datafilepath,'FreqFilter',20); 
numtrial = 1; 
data = data{numtrial}; 

%%
title_fig = "S03_Trial001_loading"; 
fig = newfig(title_fig); 
set(fig,'Position',[18.1769   13.8113   16.0338    6.6146]); 
ax = subplot_ax(1,1); hold(ax,'on');
fontsize = 20; 
set(gca, 'FontName', 'SansSerif')
frameidx = load(fullfile('analysis',"S03",'frameidx.mat')); frameidx = frameidx.frameidx; 
phase = "loading";
idx_frame_start = frameidx{2}.(phase).start(numtrial);
idx_frame_end = frameidx{2}.(phase).end(numtrial);
ax_idx = 1; 
timerange = 6.5*1e3+[idx_frame_start*1e3/50:idx_frame_end*1e3/50]; 
time = data.time; 


timescalebar = 500; %ms
y0 = 0; 
text(ax(ax_idx),time(timerange(end))-(timescalebar/1e3)/2,y0-y0*0.1,num2str(timescalebar)+"ms",...
    'FontSize',fontsize,...
    'Rotation',0,...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
plot(ax(ax_idx),[time(timerange(end))-(timescalebar/1e3) time(timerange(end))],[y0 y0],'k-','linew',2);

%Y-AXIS
x0 = time(timerange(1))-0.05;%time(timerange(end))-1.05*(time(timerange(end))-time(timerange(1))); 
xshift = 0.01;%0.05*(time(timerange(end))-time(timerange(1)));

%Y-axis - ax(1)
ax_idx = 1;
ylimplot_1 = [0,5]; 
text(ax(ax_idx),x0-10*xshift,(ylimplot_1(2)+ylimplot_1(1))/2,"Force (N)",...
        'FontSize',fontsize,...
        'Rotation',90,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom'); 

plot(ax(ax_idx),[x0 x0],ylimplot_1,'k-','linew',1.5);
yticklabelvec = ylimplot_1(1):1:ylimplot_1(2);
ytickvec = ylimplot_1(1):1:ylimplot_1(2);
for ii = 1:length(ytickvec)
    plot(ax(ax_idx),[time(timerange(1)) time(timerange(end))],[ytickvec(ii) ytickvec(ii)],'-','Color',[1 1 1].*0.5,'linew',0.5);
    plot(ax(ax_idx),[x0-3*xshift x0],[ytickvec(ii) ytickvec(ii)],'k-','linew',1.5);
    text(ax(ax_idx),x0-3*xshift,ytickvec(ii),num2str(yticklabelvec(ii)),...
        'FontSize',fontsize,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle'); 
end
set(ax(ax_idx),'Color','None'); 
idframe = [1,2,3,5,8,11,18,23]-1; 
idxframe = 6.5*1e3+idx_frame_start*1e3/50+idframe*1e3/25; 
% for ii = 1:length(idxframe)
%     text(ax(ax_idx),data.time(idxframe(ii)),data.nf(idxframe(ii)),num2str((ii)),...
%         'FontSize',15,...
%         'HorizontalAlignment','right',...
%         'VerticalAlignment','bottom');
% end
plot(ax(ax_idx),data.time(timerange), data.nf(timerange),'k-','LineWidth',5); 
plot(ax(ax_idx),data.time(idxframe), data.nf(idxframe),'s',...
    'Color','r',...
    'MarkerSize',8,...
    'MarkerFaceColor','r');%[1 1 1].*0.5);
ax(ax_idx).XAxis.Visible = 'off'; 
ax(ax_idx).YAxis.Visible = 'off'; 
%X-axis
% ax(ax_idx).XLim = [time(timerange(1)) time(timerange(end))];

