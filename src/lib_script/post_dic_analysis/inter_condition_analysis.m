%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
trial_target_set = [];
% Frictional conditions
FRcond{1} = "glass"; 
FRcond{2} = "coating";
FS_robot = 1e3; 
actual_FS_vid = 50;

%
subject = 'S09';
phase = 'loading';
material_nbr = 2;  
%
% nfcond_set = 1; 
% spddxlcond_set = 0.04; 
% r = search_trial2target(baseDataPath,subject,phase,material_nbr,...
%     nfcond_set,spddxlcond_set); 
% trial_target_set = [trial_target_set,r(end-3:end)];
%
% nfcond_set = 1; 
% spddxlcond_set = 0.08; 
% r = search_trial2target(baseDataPath,subject,phase,material_nbr,...
%     nfcond_set,spddxlcond_set); 
% trial_target_set = [trial_target_set,r(end-3:end)];
%
nfcond_set = 5; 
spddxlcond_set = 0.04; 
r = search_trial2target(baseDataPath,subject,phase,material_nbr,...
    nfcond_set,spddxlcond_set); 
trial_target_set = [trial_target_set,r(end-3:end)];
%
nfcond_set = 5; 
spddxlcond_set = 0.08; 
r = search_trial2target(baseDataPath,subject,phase,material_nbr,...
    nfcond_set,spddxlcond_set); 
trial_target_set = [trial_target_set,r(end-3:end)];

fprintf('trial target : [%s]\n',num2str(trial_target_set));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Loading trials for subject %s...\n",subject); 

Ntrial = length(trial_target_set); 

% Init outputs 
data_dic = cell(Ntrial,1); 
data_robot = cell(Ntrial,1); 

mm = 0;  
T = cell(Ntrial,1);
for itrial = trial_target_set
    mm = mm+1; 
    fprintf("--> loading %d...",itrial);
    material = FRcond{material_nbr}; 
    %set filename
    file = fullfile(baseResultPath,'analysis',subject,material,sprintf("%03d",itrial),...
        sprintf('%s_%03d_proccessed_data_%s.mat',subject,itrial,phase));
    r = load(file); r = r.data; 
    T{mm} = r;
    fprintf("done\n"); 
end
fprintf("done\n"); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RADIAL box plot fig : with FORCE condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create array for boxplot 
Nrep = 4; 
Ncond = 2; 
Nrad = length(T{1}.name_rad);
Ncirc = length(T{1}.name_circ);
val = zeros(Nrep,Ncond*Nrad,Ncirc);
faceMeasure = fieldnames(T{1}.FC);
faceMeasure = faceMeasure(1:end-1);
Nface = length(faceMeasure);

% settings : colormap, legend name,...
C=[0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
color_set_array = C(1:Nrad,:);

val_all = cell(Nface,1);
for iface = 1:Nface
faceMeasure_ii = faceMeasure{iface}; 
val = zeros(Nrep,Ncond*Nrad,Ncirc);
icount_1 = 0; 
icount_2 = 0; 
for mm = 1:Ntrial
    mask_cond = T{mm}.nfcond==1;
    val_cur = squeeze(T{mm}.FC_cum_val(end,iface,:,:));
    if mask_cond
        icount_1 = icount_1+1; 
        val(icount_1,1:Nrad,:) = val_cur;
    else 
        icount_2 = icount_2+1; 
        val(icount_2,Nrad+1:2*Nrad,:) = val_cur;
    end
end
val_all{iface} = val; 
end
cond_name = {"1N","5N"};
lgd_name = T{1}.name_rad;
ylim_set = [-15,15];

%--------------------------------------------------------------------------
% Segmentation RADIAL 
%--------------------------------------------------------------------------
% figure settings
title_fig = sprintf("%s : Radial condition NF%s-%s-%s",subject,'speckles',...
    material,phase); 
hf = newfig(title_fig);
set(hf,'Unit','Normalized','Position',[0,0.03,0.32,0.8926]);
set(hf,'Color','White') ;

% ax size
Nrow = Ncirc; 
Ncol = Nface; % boxplot total / map of regions
mergeElement={}; 
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 
for iface = 1:Nface
for xx = 1:Ncirc 
    ax_idx = iface+(xx-1)*Ncol;
    FCcurrent = squeeze(val_all{iface}(:,:,xx));
    Ndata = size(FCcurrent,1); 
    x = mat2cell(FCcurrent,Ndata,repelem(Nrad,Ncond)); 
    h = boxplotGroup(h_ax(ax_idx),x,'groupLines',true,...
        'PrimaryLabels',cond_name, ...
        'SecondaryLabels',lgd_name, ...
        'GroupLabelType','Vertical',...
        'PlotStyle','Compact','BoxStyle','filled',...
        'Colors',color_set_array,'GroupType','withinGroups');
    h_line_link = gobjects(1,Nface);
    for jj=1:Ncond
        medOuter = findobj(h.boxplotGroup(Nface+1-jj),'tag','MedianOuter');
        xpos=flip([medOuter(:).XData]);
        h_line_link(jj)=plot(h_ax(ax_idx),xpos,median(x{jj},'omitnan'),'-','Color',[1 1 1].*0.7);
    end
    %Change to color of the median points. 
    set(h.xline,'LabelVerticalAlignment','bottom')
    medOuter = findobj(h.boxplotGroup,'tag','MedianOuter');
    set(medOuter,'MarkerFaceColor','k')
    medInner = findobj(h.boxplotGroup,'tag','MedianInner');
    set(medInner,'MarkerEdgeColor','w');
    %Change outlier symbol.
    outliers = findobj(h.boxplotGroup,'Tag','Outliers');
    set(outliers, 'Marker','.'); 
    %Change width of whiskers.
    set(findobj(h.boxplotGroup,'tag','Whisker'),'LineWidth',1)
    %Rotate xtick labels and font.
    h.axis.XTickLabelRotation = 45; 
%     h.axis.XAxis.FontName = 'fixedwidth';
    h.axis.XAxis.FontWeight = 'bold';
    h.axis.XAxis.FontSize = 10;
    if xx~=Ncirc; set(h.axis.XAxis,'Visible','off'); end
    %Set the background color for each group.
    groupX = [cell2mat(get(h.xline,'value'));max(xlim(h.axis))];
    %Set YLim 
    h_ax(ax_idx).YLim = ylim_set; 
    yl = ylim(h.axis);
    ph = gobjects(1,numel(groupX)-1);
    for gg = 1:numel(groupX)-1
        ph(gg) = patch(h.axis, ...
            groupX(gg+[0,1,1,0]),...
            yl([1,1,2,2]),...
            color_set_array(gg,:),...
            'FaceAlpha',0.2,...
            'LineStyle','none');
    end
    % z-line
    h_zline = yline(h_ax(ax_idx),0,'-','LineWidth',1);
    
%     uistack(h_line_link,'bottom');
    uistack(h_zline,'bottom');
    uistack(ph,'bottom');
    %Turn on horizontal grid
    h.axis.YGrid = 'on'; 
    %Box off 
    box(h_ax(ax_idx),'off');
    h_ax(ax_idx).YLabel.FontSize = 15; 
    h_ax(ax_idx).YLabel.FontWeight = 'bold'; 
    
    if xx == 1 ; text(h_ax(ax_idx),mean(xpos)-0.5,ylim_set(2),sprintf('%s',faceMeasure{iface}),...
            'HorizontalAlignment','center','VerticalAlignment','bottom',...
            'FontSize',12,'FontWeight','bold'); 
    end
%     set(h_ax(ax_idx),'Color','None','visible','off');
%     set(findall(h_ax(ax_idx), 'type', 'text'),'visible','on')
end 
end
% allow out of plot data
set(h_ax,'Clipping','off');
%%

%--------------------------------------------------------------------------
% Segmentation CIRCUMPHERENCIAL 
%--------------------------------------------------------------------------
C = colororder(autumn(Ncirc)); %parula
color_set_array = C(1:Ncirc,:);
lgd_name = T{1}.name_circ;
color_cond{1} = [1,1,1].*0.5;
color_cond{2} = 'k';
ylim_set = [-10,10];

% figure settings
title_fig = sprintf("%s : Circumferencial condition NF%s-%s-%s",subject,'speckles',...
    material,phase); 
hf = newfig(title_fig);
set(hf,'Unit','Normalized','Position',[0,0.03,0.32,0.8926]);
set(hf,'Color','White') ;

% ax size
Nrow = Nrad; 
Ncol = Nface; % boxplot total / map of regions
mergeElement={}; 
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 
for iface = 1:Nface
for xx = 1:Nrad 
    ax_idx = iface+(xx-1)*Ncol;
    FCcurrent = val_all{iface}(:,xx+(0:Nrad:(Nrad*Ncond-1)),:);
    FCcurrent = reshape(permute(FCcurrent,[1 3 2]),Ndata,Ncond*Ncirc); 
    Ndata = size(FCcurrent,1); 
    x = mat2cell(FCcurrent,Ndata,repelem(Ncirc,Ncond)); 
   
    h = boxplotGroup(h_ax(ax_idx),x,'groupLines',true,...
        'PrimaryLabels',cond_name, ...
        'SecondaryLabels',lgd_name, ...
        'GroupLabelType','Vertical',...
        'PlotStyle','Compact','BoxStyle','filled',...
        'Colors',['k';'k';'k';'k';'k'],'GroupType','withinGroups');
    h_line_link = gobjects(1,Nface);
    for jj=1:Ncond
        medOuter = findobj(h.boxplotGroup(Nface+1-jj),'tag','MedianOuter');
        xpos=flip([medOuter(:).XData]);
        h_line_link(jj)=plot(h_ax(ax_idx),xpos,median(x{jj},'omitnan'),'-','Color',[1 1 1].*0.7);
        % overlay the scatter plots
        scatter(h_ax(ax_idx),repmat(xpos,Ndata,1),x{jj},20,"filled",'jitter','on','JitterAmount',0.2,'MarkerFaceColor',color_cond{jj});
    end
    %Change to color of the median points. 
    set(h.xline,'LabelVerticalAlignment','bottom')
    medOuter = findobj(h.boxplotGroup,'tag','MedianOuter');
    set(medOuter,'MarkerFaceColor','k')
    medInner = findobj(h.boxplotGroup,'tag','MedianInner');
    set(medInner,'MarkerEdgeColor','w');
    %Change outlier symbol.
    outliers = findobj(h.boxplotGroup,'Tag','Outliers');
    set(outliers, 'Marker','.'); 
    %Change width of whiskers.
    set(findobj(h.boxplotGroup,'tag','Whisker'),'LineWidth',1)
    %Rotate xtick labels and font.
    h.axis.XTickLabelRotation = 45; 
%     h.axis.XAxis.FontName = 'fixedwidth';
    h.axis.XAxis.FontWeight = 'bold';
    h.axis.XAxis.FontSize = 10;
    if xx~=Nrad; set(h.axis.XAxis,'Visible','off'); end
    %Set the background color for each group.
    groupX = [cell2mat(get(h.xline,'value'));max(xlim(h.axis))];
    %Set YLim 
    h_ax(ax_idx).YLim = ylim_set; 
    yl = ylim(h.axis);
    ph = gobjects(1,numel(groupX)-1);
    for gg = 1:numel(groupX)-1
        ph(gg) = patch(h.axis, ...
            groupX(gg+[0,1,1,0]),...
            yl([1,1,2,2]),...
            color_set_array(gg,:),...
            'FaceAlpha',0.2,...
            'LineStyle','none');
    end
    % z-line
    h_zline = yline(h_ax(ax_idx),0,'-','LineWidth',1);
    
%     uistack(h_line_link,'bottom');
    uistack(h_zline,'bottom');
    uistack(ph,'bottom');
    %Turn on horizontal grid
    h.axis.YGrid = 'on'; 
    %Box off 
    box(h_ax(ax_idx),'off');
    h_ax(ax_idx).YLabel.FontSize = 15; 
    h_ax(ax_idx).YLabel.FontWeight = 'bold'; 
    
    if xx == 1 ; text(h_ax(ax_idx),mean(xpos)-0.5,ylim_set(2),sprintf('%s',faceMeasure{iface}),...
            'HorizontalAlignment','center','VerticalAlignment','bottom',...
            'FontSize',12,'FontWeight','bold'); 
    end
%     set(h_ax(ax_idx),'Color','None','visible','off');
%     set(findall(h_ax(ax_idx), 'type', 'text'),'visible','on')
end 
end
% allow out of plot data
set(h_ax,'Clipping','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RADIAL box plot fig : with SPD condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create array for boxplot 
Nrep = 4; 
Ncond = 2; 
Nrad = length(T{1}.name_rad);
Ncirc = length(T{1}.name_circ);
val = zeros(Nrep,Ncond*Nrad,Ncirc);
faceMeasure = fieldnames(T{1}.FC);
faceMeasure = faceMeasure(1:end-1);
Nface = length(faceMeasure);

% settings : colormap, legend name,...
C=[0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...	
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
color_set_array = C(1:Nrad,:);

val_all = cell(Nface,1);
for iface = 1:Nface
faceMeasure_ii = faceMeasure{iface}; 
val = zeros(Nrep,Ncond*Nrad,Ncirc);
icount_1 = 0; 
icount_2 = 0; 
for mm = 1:Ntrial
    mask_cond = T{mm}.spddxlcond==0.04;
    val_cur = squeeze(T{mm}.FC_cum_val(end,iface,:,:));
    if mask_cond
        icount_1 = icount_1+1; 
        val(icount_1,1:Nrad,:) = val_cur;
    else 
        icount_2 = icount_2+1; 
        val(icount_2,Nrad+1:2*Nrad,:) = val_cur;
    end
end
val_all{iface} = val; 
end

cond_name = {"low","high"};

%--------------------------------------------------------------------------
% Segmentation CIRCUMPHERENCIAL 
%--------------------------------------------------------------------------
C = colororder(autumn(Ncirc)); %parula
color_set_array = C(1:Ncirc,:);
lgd_name = T{1}.name_circ;
color_cond{1} = [1,1,1].*0.5;
color_cond{2} = 'k';
ylim_set = [-10,10];

% figure settings
title_fig = sprintf("%s : Circumferencial condition spd%s-%s-%s",subject,'speckles',...
    material,phase); 
hf = newfig(title_fig);
set(hf,'Unit','Normalized','Position',[0,0.03,0.32,0.8926]);
set(hf,'Color','White') ;

% ax size
Nrow = Nrad; 
Ncol = Nface; % boxplot total / map of regions
mergeElement={}; 
h_ax = subplot_ax(Nrow,Ncol,'merge',mergeElement,'shownum',0); hold(h_ax,'on'); 
for iface = 1:Nface
for xx = 1:Nrad 
    ax_idx = iface+(xx-1)*Ncol;
    FCcurrent = val_all{iface}(:,xx+(0:Nrad:(Nrad*Ncond-1)),:);
    Ndata = size(FCcurrent,1); 
    FCcurrent = reshape(permute(FCcurrent,[1 3 2]),Ndata,Ncond*Ncirc); 
    x = mat2cell(FCcurrent,Ndata,repelem(Ncirc,Ncond)); 
   
    h = boxplotGroup(h_ax(ax_idx),x,'groupLines',true,...
        'PrimaryLabels',cond_name, ...
        'SecondaryLabels',lgd_name, ...
        'GroupLabelType','Vertical',...
        'PlotStyle','Compact','BoxStyle','filled',...
        'Colors',['k';'k';'k';'k';'k'],'GroupType','withinGroups');
    h_line_link = gobjects(1,Nface);
    for jj=1:Ncond
        medOuter = findobj(h.boxplotGroup(Nface+1-jj),'tag','MedianOuter');
        xpos=flip([medOuter(:).XData]);
        h_line_link(jj)=plot(h_ax(ax_idx),xpos,median(x{jj},'omitnan'),'-','Color',[1 1 1].*0.7);
        % overlay the scatter plots
        scatter(h_ax(ax_idx),repmat(xpos,Ndata,1),x{jj},20,"filled",'jitter','on','JitterAmount',0.2,'MarkerFaceColor',color_cond{jj});
    end
    %Change to color of the median points. 
    set(h.xline,'LabelVerticalAlignment','bottom')
    medOuter = findobj(h.boxplotGroup,'tag','MedianOuter');
    set(medOuter,'MarkerFaceColor','k')
    medInner = findobj(h.boxplotGroup,'tag','MedianInner');
    set(medInner,'MarkerEdgeColor','w');
    %Change outlier symbol.
    outliers = findobj(h.boxplotGroup,'Tag','Outliers');
    set(outliers, 'Marker','.'); 
    %Change width of whiskers.
    set(findobj(h.boxplotGroup,'tag','Whisker'),'LineWidth',1)
    %Rotate xtick labels and font.
    h.axis.XTickLabelRotation = 45; 
%     h.axis.XAxis.FontName = 'fixedwidth';
    h.axis.XAxis.FontWeight = 'bold';
    h.axis.XAxis.FontSize = 10;
    if xx~=Nrad; set(h.axis.XAxis,'Visible','off'); end
    %Set the background color for each group.
    groupX = [cell2mat(get(h.xline,'value'));max(xlim(h.axis))];
    %Set YLim 
    h_ax(ax_idx).YLim = ylim_set; 
    yl = ylim(h.axis);
    ph = gobjects(1,numel(groupX)-1);
    for gg = 1:numel(groupX)-1
        ph(gg) = patch(h.axis, ...
            groupX(gg+[0,1,1,0]),...
            yl([1,1,2,2]),...
            color_set_array(gg,:),...
            'FaceAlpha',0.2,...
            'LineStyle','none');
    end
    % z-line
    h_zline = yline(h_ax(ax_idx),0,'-','LineWidth',1);
    
%     uistack(h_line_link,'bottom');
    uistack(h_zline,'bottom');
    uistack(ph,'bottom');
    %Turn on horizontal grid
    h.axis.YGrid = 'on'; 
    %Box off 
    box(h_ax(ax_idx),'off');
    h_ax(ax_idx).YLabel.FontSize = 15; 
    h_ax(ax_idx).YLabel.FontWeight = 'bold'; 
    
    if xx == 1 ; text(h_ax(ax_idx),mean(xpos)-0.5,ylim_set(2),sprintf('%s',faceMeasure{iface}),...
            'HorizontalAlignment','center','VerticalAlignment','bottom',...
            'FontSize',12,'FontWeight','bold'); 
    end
%     set(h_ax(ax_idx),'Color','None','visible','off');
%     set(findall(h_ax(ax_idx), 'type', 'text'),'visible','on')
end 
end
% allow out of plot data
set(h_ax,'Clipping','off');

