 
baseDataPath = "E:\transfert_SMP"; 
datafilepath = fullfile(baseDataPath,"rawdata","S03","speckles","coating","robot");
data = import_robotdata_subject(datafilepath,'FreqFilter',20); 
numtrial = 9; 
data = data{numtrial}; 

%%
title_fig = "S01_Trial009_tangential"; 
fig = newfig(title_fig); 
set(fig,'Position',[18.1769   13.8113   16.0338    6.6146]); 
ax = subplot_ax(1,1); hold(ax,'on');
fontsize = 20; 
set(gca, 'FontName', 'SansSerif')
frameidx = load(fullfile('analysis',"S03",'frameidx.mat')); frameidx = frameidx.frameidx; 
phase = "loading";
ax_idx = 1; 
timerange = (6.5+3)*1e3:(6.5+3+1.6)*1e3; 
time = data.time; 

y0 = -0.1; 
yshift = 0.1;  
xlimplot_1 = [0,8]; 
text(ax(ax_idx),(xlimplot_1(2)+xlimplot_1(1))/2,y0-10*yshift,"Position (mm)",...
        'FontSize',fontsize,...
        'Rotation',0,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top'); 
plot(ax(ax_idx),xlimplot_1,[y0, y0],'k-','linew',1.5);
xticklabelvec = xlimplot_1(1):2:xlimplot_1(2);
xtickvec = xlimplot_1(1):2:xlimplot_1(2);
for ii = 1:length(xtickvec)
    plot(ax(ax_idx),[xtickvec(ii) xtickvec(ii)],[y0-3*yshift y0],'k-','linew',1.5);
    text(ax(ax_idx),xtickvec(ii),y0-3*yshift,num2str(xticklabelvec(ii)),...
        'FontSize',fontsize,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top'); 
end

%Y-AXIS
x0 = -0.3; %time(timerange(1))-0.05;%time(timerange(end))-1.05*(time(timerange(end))-time(timerange(1))); 
xshift = 0.05;%0.05*(time(timerange(end))-time(timerange(1)));

%Y-axis - ax(1)
ax_idx = 1;
ylimplot_1 = [0,6]; 
text(ax(ax_idx),x0-10*xshift,(ylimplot_1(2)+ylimplot_1(1))/2,"Force (N)",...
        'FontSize',fontsize,...
        'Rotation',90,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom'); 

plot(ax(ax_idx),[x0 x0],ylimplot_1,'k-','linew',1.5);
yticklabelvec = ylimplot_1(1):2:ylimplot_1(2);
ytickvec = ylimplot_1(1):2:ylimplot_1(2);
for ii = 1:length(ytickvec)
    if ii ~=1 
        plot(ax(ax_idx),xlimplot_1,[ytickvec(ii) ytickvec(ii)],'-','Color',[1 1 1].*0.5,'linew',0.5);
    end
   plot(ax(ax_idx),[x0-3*xshift x0],[ytickvec(ii) ytickvec(ii)],'k-','linew',1.5);
    text(ax(ax_idx),x0-3*xshift,ytickvec(ii),num2str(yticklabelvec(ii)),...
        'FontSize',fontsize,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle'); 
end
set(ax(ax_idx),'Color','None'); 
idframe = [1,8,12,16,20,26,31,35,38]-1; 
idxframe = (6.5+3)*1e3+idframe*1e3/25; 
% for ii = 1:length(idxframe)
%     text(ax(ax_idx),data.time(idxframe(ii)),data.nf(idxframe(ii)),num2str((ii)),...
%         'FontSize',15,...
%         'HorizontalAlignment','right',...
%         'VerticalAlignment','bottom');
% end
pos = abs(data.position_robot(:,1)-data.position_robot(timerange(1),1)); 
plot(ax(ax_idx),[pos(timerange(1)),pos(timerange(end))], [5 5] ,'--','LineWidth',1,'Color',[1 1 1].*0.5); 
plot(ax(ax_idx),pos(timerange), data.nf(timerange),'-','LineWidth',5,'Color',[1 1 1].*0.5); 
plot(ax(ax_idx),pos(timerange), data.tf(timerange),'k-','LineWidth',5); 
plot(ax(ax_idx),pos(idxframe), data.tf(idxframe),'s',...
    'Color','r',...
    'MarkerSize',8,...
    'MarkerFaceColor','r');%[1 1 1].*0.5);
ax(ax_idx).XAxis.Visible = 'off'; 
ax(ax_idx).YAxis.Visible = 'off'; 
%X-axis
% ax(ax_idx).XLim = [time(timerange(1)) time(timerange(end))];

%legend
qw = cell(2,1); 
qw{1} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],[1 1 1].*0);
qw{2} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],[1 1 1].*0.5);
hlegend = legend(ax(ax_idx),[qw{:}],{'TF','NF'}); 
set(hlegend, 'NumColumns', 2,...
    'Location','southeast',...
    'Box','off',...
    'FontSize',15);
box(hlegend,'off'); 


