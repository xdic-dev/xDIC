

subject_set = ["S05"];%,"S06","S07","S08","S09"];
table_normal_loading = table(); 
for ii = 1:length(subject_set)
    temptbl = load(fullfile('analysis',sprintf('%s_table_normal_loading',subject_set(ii))));
    table_normal_loading = [table_normal_loading;temptbl.table_normal_loading];
end
fprintf('done\n'); 
table_normal_loading(1:10,:)

subject_list = unique(table_normal_loading.Subject);
spdcond_list = unique(table_normal_loading.SPDcond); 
nfcond_list = unique(table_normal_loading.NF_formatted); 

% Correction values
table_normal_loading.Area(isnan(table_normal_loading.Area)) = 0; 
for isubject = subject_list'
    mask_target = strcmp(table_normal_loading.Subject,isubject);
    if str2double(cell2mat(regexp(isubject, '\d+', 'match')))<8
        table_normal_loading.SPDcond(mask_target & (table_normal_loading.TrialNbr<=24)) = 0.04; 
        table_normal_loading.SPDcond(mask_target & (table_normal_loading.TrialNbr>24)) = 0.08; 
    end
end

%% Trial averaged graphs: 

dummyfun1 = @(x)median(x,'omitnan'); 
dummyfun2 = @(x)min(x,[],'omitnan'); 
dummyfun3 = @(x)max(x,[],'omitnan'); 
dummyfun4 = @(x)mean(x,'omitnan'); 
dummyfun5 = @(x)std(x,'omitnan'); 
dummyfun6 = @(x)sum(x,'omitnan'); 

%Define the strain deformation type that you want to plot and traces at the
%normal force wanted. 
def_type = "cum"; %rate %cum 
nf_target = [0.1,0.2,0.5,1,5];%[0.2,1,5];%nfcond_list(end); 
color_spd = {'b','r'}; 

%--------------------------------------------------------------------------
% CONTACT AREA vs subject vs speed
groupvars = ["Subject","SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5,dummyfun6}; 
targetvar = ["Area"]; 
tbl = table_normal_loading;
tblstats0 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats0(1:5,:)


ylim_set = [0,250];
title_fig = sprintf("Contact_area"); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.0063    0.6181    0.2547    0.2162]); 
Nrow = 1;length(spdcond_list); Ncol = length(subject_list);
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',10,'Clipping','off','Color','none'); 
for ii = 1:length(subject_list)
subject_ii = subject_list(ii); 
for jj = 1:length(spdcond_list)
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.SegFactor == tbl.SegFactor(end)...
        & tbl.SPDcond == spdcond_jj...
        & tbl.Subject == subject_ii; 
    logic_in_cur = tbl.logic_in(mask_target);
    nfcond_cur = tbl.NFcond(mask_target);
    x1 = tbl.NF_formatted(mask_target);
    y1 = tbl.(targetvar(1))(mask_target);
    
    mask_target = tblstats0.SegFactor == tblstats0.SegFactor(end)...
        & tblstats0.SPDcond == spdcond_jj...
        & tblstats0.Subject == subject_ii; 


    x2 = tblstats0.NF_formatted(mask_target);
    y3 = tblstats0.('Fun1_'+targetvar(1))(mask_target); 
    y5 = tblstats0.('Fun5_'+targetvar(1))(mask_target); 
 
    %graph 
    ax_idx = ii;%+(jj-1)*Ncol; 
    h1 = yline(h_ax(ax_idx),0);
%     scatter(h_ax(ax_idx),x1,y1,20,[1 1 1].*0.5,"filled"); 

    h_data = plot(h_ax(ax_idx), x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    xline(h_ax(ax_idx),nf_target,'Color',[1 1 1].*0.5); 
    mask_target = logical(sum(x2==nf_target,2));
    hpoint = plot(h_ax(ax_idx),x2(mask_target),y3(mask_target),...
        'Marker','.','LineStyle','none','MarkerSize',20,'Color',[1 1 1].*0.5);
    
    if ii == 1 
        ylabel(h_ax(ax_idx),'Area (mm2)');
    else
        h_ax(ax_idx).YTickLabel = [];
    end
    xlabel(h_ax(ax_idx),'Force (N)');
    
    xlim(h_ax(ax_idx),[0 nfcond_cur(1)]); 
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim'); ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == 1
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("%s\nspd = %1.2f",subject_ii,spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
    uistack(h2,'bottom');
    uistack(h1,'bottom');
%     uistack(h4,'bottom');
end
end
qw = cell(1,length(spdcond_list));
for jj = 1:length(spdcond_list)
qw{jj} = plot(h_ax(Ncol),NaN,NaN,'-','Color',color_spd{jj},'LineWidth',2);
end
h_leg = legend(h_ax(Ncol),[hpoint,qw{:}],{'eval point','low speed','high speed'});
set(h_leg,'Position',[0.3940    0.9359    0.2311    0.0493],'box','on',...
    'NumColumns',3,'Color','white'); 

%
%--------------------------------------------------------------------------
% CONTACT AREA NORMALIZED vs speed
groupvars = ["SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5,dummyfun6}; 
targetvar = ["Area"]; 
%normalized area 
temptbl = grpstats(table_normal_loading,"Subject",whichstats,"DataVars","Area"); 
tbl = table_normal_loading;
for isubject = subject_list'
    mask_target = strcmp(tbl.Subject,isubject); 
    tbl.Area(mask_target) = tbl.Area(mask_target)/temptbl.Fun3_Area(strcmp(temptbl.Subject,isubject));
end
tblstats0 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats0(1:5,:)


ylim_set = [0,1];
title_fig = sprintf("Contact_area_all"); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.2615    0.6185    0.1271    0.2162]); 
Nrow = 1;length(spdcond_list); Ncol = 1;
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',10,'Clipping','off'); 
for jj = 1:length(spdcond_list)
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.SegFactor == tbl.SegFactor(1)...
        & tbl.SPDcond == spdcond_jj; 
    logic_in_cur = tbl.logic_in(mask_target);
    
    x1 = tbl.NF_formatted(mask_target);
    y1 = tbl.(targetvar(1))(mask_target); 
    
    mask_target = tblstats0.SegFactor == tblstats0.SegFactor(1)...
        & tblstats0.SPDcond == spdcond_jj; 


    x2 = tblstats0.NF_formatted(mask_target);
    y3 = tblstats0.('Fun1_'+targetvar(1))(mask_target); 
    y5 = tblstats0.('Fun5_'+targetvar(1))(mask_target); 
 
    %graph 
    ax_idx = 1;jj; 
    h1 = yline(h_ax(ax_idx),0);
%     scatter(h_ax(ax_idx),x1,y1,20,[1 1 1].*0.5,"filled"); 

    plot(h_ax(ax_idx), x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    for kk = 1:length(nf_target)
        xline(h_ax(ax_idx),nf_target(kk),'Color',[1 1 1].*0.5); 
        plot(h_ax(ax_idx),x2(x2==nf_target(kk)),y3(x2==nf_target(kk)),'.','MarkerSize',20,'Color',[1 1 1].*0.5);
    end
     
    ylabel(h_ax(ax_idx),'Norm. area (-)');
    xlabel(h_ax(ax_idx),'Force (N)');
    
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == 1
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("%s\nspd = %1.2f",subject_ii,spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
    uistack(h2,'bottom');
    uistack(h1,'bottom'); 
%     uistack(h4,'bottom');
end

%%
%--------------------------------------------------------------------------
% TANGENTIAL FORCE vs subject vs speed
groupvars = ["Subject","SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5,dummyfun6}; 
targetvar = ["force_ratio"]; 
tbl = table_normal_loading;
tbl.force_ratio = tbl.TF_raw./tbl.NF_raw; 
tblstats0 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats0(1:5,:)

ylim_set = [0,3];
title_fig = sprintf("Friction estimation"); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.0063    0.6181    0.2547    0.2162]); 
Nrow = 1;length(spdcond_list); Ncol = length(subject_list);
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',10,'Clipping','off','Color','none'); 
for ii = 1:length(subject_list)
subject_ii = subject_list(ii); 
for jj = 1:length(spdcond_list)
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.SegFactor == tbl.SegFactor(end)...
        & tbl.SPDcond == spdcond_jj...
        & tbl.Subject == subject_ii; 
    logic_in_cur = tbl.logic_in(mask_target);
    trial_cur = tbl.TrialNbr(mask_target);
    nfcond_cur = tbl.NFcond(mask_target);
    x1 = tbl.NF_formatted(mask_target);
    y1 = tbl.(targetvar(1))(mask_target);
    
    mask_target = tblstats0.SegFactor == tblstats0.SegFactor(end)...
        & tblstats0.SPDcond == spdcond_jj...
        & tblstats0.Subject == subject_ii; 


    x2 = tblstats0.NF_formatted(mask_target);
    y3 = tblstats0.('Fun1_'+targetvar(1))(mask_target); 
    y5 = tblstats0.('Fun5_'+targetvar(1))(mask_target); 
 
    %graph 
    ax_idx = ii;%+(jj-1)*Ncol; 
    h1 = yline(h_ax(ax_idx),0);
    for itrial = unique(trial_cur)'
        mask_target = trial_cur == itrial; 
        plot(h_ax(ax_idx),x1(mask_target),y1(mask_target),'-','Color',color_spd{jj}); 
    end
    scatter(h_ax(ax_idx),x1,y1,10,[1 1 1].*0.5,"filled"); 

    h_data = plot(h_ax(ax_idx), x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    
    if ii == 1 
        ylabel(h_ax(ax_idx),'TF/NF (-)');
    else
        h_ax(ax_idx).YTickLabel = [];
    end
    xlabel(h_ax(ax_idx),'Force (N)');
    
    xlim(h_ax(ax_idx),[0 nfcond_cur(1)]); 
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim'); ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == 1
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("%s\nspd = %1.2f",subject_ii,spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
    uistack(h2,'bottom');
    uistack(h1,'bottom');
%     uistack(h4,'bottom');
end
end
qw = cell(1,length(spdcond_list));
for jj = 1:length(spdcond_list)
qw{jj} = plot(h_ax(Ncol),NaN,NaN,'-','Color',color_spd{jj},'LineWidth',2);
end
h_leg = legend(h_ax(Ncol),[qw{:}],{'low speed','high speed'});
set(h_leg,'Position',[0.3940    0.9359    0.2311    0.0493],'box','on',...
    'NumColumns',3,'Color','white'); 

%%
%--------------------------------------------------------------------------
% CONTACT AREA NORMALIZED vs speed
groupvars = ["SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5,dummyfun6}; 
targetvar = ["Area"]; 
%normalized area 
temptbl = grpstats(table_normal_loading,"Subject",whichstats,"DataVars","Area"); 
tbl = table_normal_loading;
for isubject = subject_list'
    mask_target = strcmp(tbl.Subject,isubject); 
    tbl.Area(mask_target) = tbl.Area(mask_target)/temptbl.Fun3_Area(strcmp(temptbl.Subject,isubject));
end
tblstats0 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats0(1:5,:)


ylim_set = [0,1];
title_fig = sprintf("Contact_area_all"); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.2615    0.6185    0.1271    0.2162]); 
Nrow = 1;length(spdcond_list); Ncol = 1;
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',10,'Clipping','off'); 
for jj = 1:length(spdcond_list)
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.SegFactor == tbl.SegFactor(1)...
        & tbl.SPDcond == spdcond_jj; 
    logic_in_cur = tbl.logic_in(mask_target);
    
    x1 = tbl.NF_formatted(mask_target);
    y1 = tbl.(targetvar(1))(mask_target); 
    
    mask_target = tblstats0.SegFactor == tblstats0.SegFactor(1)...
        & tblstats0.SPDcond == spdcond_jj; 


    x2 = tblstats0.NF_formatted(mask_target);
    y3 = tblstats0.('Fun1_'+targetvar(1))(mask_target); 
    y5 = tblstats0.('Fun5_'+targetvar(1))(mask_target); 
 
    %graph 
    ax_idx = 1;jj; 
    h1 = yline(h_ax(ax_idx),0);
%     scatter(h_ax(ax_idx),x1,y1,20,[1 1 1].*0.5,"filled"); 

    plot(h_ax(ax_idx), x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    for kk = 1:length(nf_target)
        xline(h_ax(ax_idx),nf_target(kk),'Color',[1 1 1].*0.5); 
        plot(h_ax(ax_idx),x2(x2==nf_target(kk)),y3(x2==nf_target(kk)),'.','MarkerSize',20,'Color',[1 1 1].*0.5);
    end
     
    ylabel(h_ax(ax_idx),'Norm. area (-)');
    xlabel(h_ax(ax_idx),'Force (N)');
    
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == 1
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("%s\nspd = %1.2f",subject_ii,spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
    uistack(h2,'bottom');
    uistack(h1,'bottom');
%     uistack(h4,'bottom');
end


%--------------------------------------------------------------------------
% PRINCIPAL STRAINS 
% - Col : subject
% - Row : speed
groupvars = ["Subject","SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5}; 
targetvar = ["Epc1_"+def_type,"Epc2_"+def_type,"logic_in"]; 
tbl = table_normal_loading;
tblstats1 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats1(1:5,:)

% graph 

if strcmp(def_type,"cum")
    title_def_type = "Total"; 
    ylim_set = [-10,10];
    ylabel_set = "Total strain (%)";
elseif strcmp(def_type,"rate")
    title_def_type = "Rate"; 
    ylim_set = [-70,70];
    ylabel_set = "Rate strain (%/s)";
end

title_fig = sprintf("%s_%1.2f",title_def_type,nf_target); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.3948    0.4019    0.2547    0.4324]); 
Nrow = length(spdcond_list); Ncol = length(subject_list);
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',12,'Clipping','off'); 
for ii = 1:length(subject_list)
subject_ii = subject_list(ii); 
for jj = 1:length(spdcond_list)
h_pin1 = gobjects(1,length(nf_target));
h_pin2 = gobjects(1,length(nf_target));
iplot = 0; 
for kk = 1:length(nf_target)
    iplot = iplot+1; 
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.NF_formatted == nf_target(kk)...
        & tbl.SPDcond == spdcond_jj...
        & tbl.Subject == subject_ii; 
    logic_in_cur = tbl.logic_in(mask_target);
    trials_number_cur = tbl.TrialNbr(mask_target);
    x1 = tbl.SegFactor(mask_target);
    y1 = tbl.(targetvar(1))(mask_target); 
    y2 = tbl.(targetvar(2))(mask_target); 

    mask_target = tblstats1.NF_formatted == nf_target(kk)...
        & tblstats1.SPDcond == spdcond_jj...
        & tblstats1.Subject == subject_ii; 

    idx_logic = tblstats1.('Fun4_'+targetvar(3))(mask_target)>=0.5;
    idx_logic_first = find(~idx_logic(2:end),1,'first');
    x2 = tblstats1.SegFactor(mask_target);
    y3 = tblstats1.('Fun4_'+targetvar(1))(mask_target); 
    y4 = tblstats1.('Fun4_'+targetvar(2))(mask_target); 
    y5 = tblstats1.('Fun5_'+targetvar(1))(mask_target); 
    y6 = tblstats1.('Fun5_'+targetvar(2))(mask_target); 
    
    ax_idx = ii+(jj-1)*Ncol; 
    h1 = yline(h_ax(ax_idx),0);
    %indiv trials
    for itrials = unique(tbl.TrialNbr)'
        plot(h_ax(ax_idx),x1(trials_number_cur==itrials),y1(trials_number_cur==itrials),'Color',[1 1 1].*0.8,'LineWidth',1); 
        plot(h_ax(ax_idx),x1(trials_number_cur==itrials),y2(trials_number_cur==itrials),'Color',[1 1 1].*0.8,'LineWidth',1); 
    end
    scatter(h_ax(ax_idx),x1(~logic_in_cur),y1(~logic_in_cur),10,[1 1 1].*0.5,'Marker','s','MarkerFaceColor','w'); 
    scatter(h_ax(ax_idx),x1(~logic_in_cur),y2(~logic_in_cur),10,[1 1 1].*0.5,'Marker','s','MarkerFaceColor','w'); 
    scatter(h_ax(ax_idx),x1(logic_in_cur),y1(logic_in_cur),10,[1 1 1].*0.5,"filled"); 
    scatter(h_ax(ax_idx),x1(logic_in_cur),y2(logic_in_cur),10,[1 1 1].*0.5,"filled"); 
    
    %averaged trials
    plot(h_ax(ax_idx), x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    plot(h_ax(ax_idx), x2,y4,'-','Color',color_spd{jj},'LineWidth',2);
    h_pin1(iplot) = plot(h_ax(ax_idx),x2(idx_logic_first),y3(idx_logic_first),'.','Color','k','MarkerSize',20); 
    h_pin2(iplot) = plot(h_ax(ax_idx),x2(idx_logic_first),y4(idx_logic_first),'.','Color','k','MarkerSize',20); 
%     h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
%     h3 = ciplot(y4-y6,y4+y6,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    
    %graph info 
    if ii == 1 
        ylabel(h_ax(ax_idx),ylabel_set);
    else
        set(get(h_ax(ax_idx), 'YAxis'), 'Visible', 'off'); 
%         h_ax(ax_idx).YTickLabel = [];
    end
    if jj ~= length(spdcond_list)
        set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off'); 
    else
        xlabel(h_ax(ax_idx),'Norm. meridional contact');
    end
    
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == 1
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("%s\nspd = %1.2f",subject_ii,spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
%     uistack(h3,'bottom');
    uistack(h2,'bottom');
    uistack(h1,'bottom');
%     uistack(h4,'bottom');
end
uistack(h_pin1,'top');
uistack(h_pin2,'top');
end
end

%--------------------------------------------------------------------------
% PRINCIPAL STRAINS - TRIAL AVERAGED : 
% - Row : speed
% Same but gathered across subjects

groupvars = ["SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5}; 
targetvar = ["Epc1_"+def_type,"Epc2_"+def_type,"logic_in"]; 
tbl = table_normal_loading;
tbl(1:5,:)
tblstats2 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats2(1:5,:)

if strcmp(def_type,"cum")
    title_def_type = "Total"; 
    ylim_set = [-10,10];
    ylabel_set = "Total strain (%)";
elseif strcmp(def_type,"rate")
    title_def_type = "Rate"; 
    ylim_set = [-70,70];
    ylabel_set = "Rate strain (%/s)";
end

title_fig = sprintf("%s_%1.2f_allsubject",title_def_type,nf_target); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.6505    0.4019    0.1271    0.4324]); 
Nrow = length(spdcond_list); Ncol = 1;
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',12,'Clipping','off'); 

for jj = 1:length(spdcond_list)
h_pin1 = gobjects(1,length(nf_target));
h_pin2 = gobjects(1,length(nf_target));
iplot = 0; 
for kk = 1:length(nf_target)
    iplot = iplot+1; 
    
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.NF_formatted == nf_target(kk)...
        & tbl.SPDcond == spdcond_jj; 

    
    x1 = tbl.SegFactor(mask_target);
    y1 = tbl.(targetvar(1))(mask_target); 
    y2 = tbl.(targetvar(2))(mask_target); 

    mask_target = tblstats2.NF_formatted == nf_target(kk)...
        & tblstats2.SPDcond == spdcond_jj; 


    idx_logic = tblstats2.('Fun4_'+targetvar(3))(mask_target)>=0.5;
    idx_logic_first = find(~idx_logic(2:end),1,'first');
    x2 = tblstats2.SegFactor(mask_target);
    y3 = tblstats2.('Fun4_'+targetvar(1))(mask_target); 
    y4 = tblstats2.('Fun4_'+targetvar(2))(mask_target); 
    y5 = tblstats2.('Fun5_'+targetvar(1))(mask_target); 
    y6 = tblstats2.('Fun5_'+targetvar(2))(mask_target); 
 
    %graph 
    ax_idx = 1+(jj-1)*Ncol; 
    h1 = yline(h_ax(ax_idx),0);
    plot(h_ax(ax_idx),x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    plot(h_ax(ax_idx),x2,y4,'-','Color',color_spd{jj},'LineWidth',2);
    h_pin1(iplot) = plot(h_ax(ax_idx),x2(idx_logic_first),y3(idx_logic_first),'.','Color','k','MarkerSize',20); 
    h_pin2(iplot) = plot(h_ax(ax_idx),x2(idx_logic_first),y4(idx_logic_first),'.','Color','k','MarkerSize',20); 
    h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    h3 = ciplot(y4-y6,y4+y6,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    
       
    ylabel(h_ax(ax_idx),ylabel_set);
    if jj ~= length(spdcond_list) && kk == length(nf_target)
        set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off'); 
    else
        xlabel(h_ax(ax_idx),'Norm. meridional contact');
    end
    
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == length(nf_target)
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("spd = %1.2f",spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
    uistack(h3,'bottom');
    uistack(h2,'bottom');
    uistack(h1,'bottom');
    
%     uistack(h4,'bottom');
end
uistack(h_pin1,'top');
uistack(h_pin2,'top');
end

title_fig = sprintf("%s_%1.2f_allsubject_spdcombined",title_def_type,nf_target); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.7781    0.6185    0.1271    0.2157]); 
Nrow = 1;length(spdcond_list); Ncol = 1;
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); grid(h_ax,'on');
set(h_ax,'FontSize',12,'Clipping','off'); 

for jj = 1:length(spdcond_list)
h_pin1 = gobjects(1,length(nf_target));
h_pin2 = gobjects(1,length(nf_target));
iplot = 0; 
for kk = 1:length(nf_target)
    iplot = iplot+1; 
    
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.NF_formatted == nf_target(kk)...
        & tbl.SPDcond == spdcond_jj; 

    
    x1 = tbl.SegFactor(mask_target);
    y1 = tbl.(targetvar(1))(mask_target); 
    y2 = tbl.(targetvar(2))(mask_target); 

    mask_target = tblstats2.NF_formatted == nf_target(kk)...
        & tblstats2.SPDcond == spdcond_jj; 


    idx_logic = tblstats2.('Fun4_'+targetvar(3))(mask_target)>=0.5;
    idx_logic_first = find(~idx_logic(2:end),1,'first');
    x2 = tblstats2.SegFactor(mask_target);
    y3 = tblstats2.('Fun4_'+targetvar(1))(mask_target); 
    y4 = tblstats2.('Fun4_'+targetvar(2))(mask_target); 
    y5 = tblstats2.('Fun5_'+targetvar(1))(mask_target); 
    y6 = tblstats2.('Fun5_'+targetvar(2))(mask_target); 
 
    %graph 
    ax_idx = 1;%+(jj-1)*Ncol; 
    h1 = yline(h_ax(ax_idx),0);
    plot(h_ax(ax_idx),x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    plot(h_ax(ax_idx),x2,y4,'-','Color',color_spd{jj},'LineWidth',2);
    h_pin1(iplot) = plot(h_ax(ax_idx),x2(idx_logic_first),y3(idx_logic_first),'.','Color','k','MarkerSize',20); 
    h_pin2(iplot) = plot(h_ax(ax_idx),x2(idx_logic_first),y4(idx_logic_first),'.','Color','k','MarkerSize',20); 
    h2 = ciplot(y3-y5,y3+y5,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    h3 = ciplot(y4-y6,y4+y6,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    
       
    ylabel(h_ax(ax_idx),ylabel_set);
    if jj ~= length(spdcond_list) && kk == length(nf_target)
        set(get(h_ax(ax_idx), 'XAxis'), 'Visible', 'off'); 
    else
        xlabel(h_ax(ax_idx),'Norm. meridional contact');
    end
    
    ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    if kk == length(nf_target)
        text(h_ax(ax_idx),mean(xlim_get),ylim_get(2),...
            sprintf("spd = %1.2f",spdcond_jj),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
    end
%     h4 = ciplot([ylim_get(1),ylim_get(1)],[ylim_get(2),ylim_get(2)],[0,1],'y',h_ax(ax_idx),0.15);
    uistack(h3,'bottom');
    uistack(h2,'bottom');
    uistack(h1,'bottom');
    
%     uistack(h4,'bottom');
end
uistack(h_pin1,'top');
uistack(h_pin2,'top');
end
%%
% graph 2 

def_type = "rate"; %rate %cum 

groupvars = ["Subject","SegFactor","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5}; 
targetvar = ["Epc1_"+def_type,"Epc2_"+def_type,"logic_in"]; 
tbl = table_normal_loading;
tblstats1 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats1(1:5,:)

groupvars = ["Subject","NF_formatted","SPDcond"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5}; 
targetvar = ["Fun4_"+targetvar(1),"Fun4_"+targetvar(2),"Fun4_"+targetvar(3)]; 
tbl = tblstats1;
tblstats2 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats2(1:20,:)

ylim_set = [-80,80];
title_fig = sprintf("MaxRate_vs_force"); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.0734    0.0287+0.8926/2    0.2641    0.8926/2]); 
Nrow = 1; Ncol = length(subject_list);
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
set(h_ax,'FontSize',12,'Clipping','off'); 
for ii = 1:length(subject_list)
subject_ii = subject_list(ii); 
for jj = 1:length(spdcond_list)
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tbl.SPDcond == spdcond_jj...
        & tbl.Subject == subject_ii; 

    x1 = tbl.NF_formatted(mask_target);
    y1 = tbl.(targetvar(1))(mask_target); 
    y2 = tbl.(targetvar(2))(mask_target); 
    idx_logic = tbl.(targetvar(3))(mask_target)>=0.5; 

    mask_target = tblstats2.SPDcond == spdcond_jj...
        & tblstats2.Subject == subject_ii; 

    x2 = tblstats2.NF_formatted(mask_target);
    y3 = tblstats2.('Fun2_'+targetvar(1))(mask_target); 
    y4 = tblstats2.('Fun3_'+targetvar(2))(mask_target); 

    mask_target = tbl.SPDcond == spdcond_jj...
        & tbl.Subject == subject_ii;
 
    %graph 
    ax_idx = ii; 
    h1 = yline(h_ax(ax_idx),0);
%     scatter(h_ax(ax_idx), x1(~idx_logic),y1(~idx_logic),20,[1 1 1].*0.5,'Marker','s','MarkerFaceColor','w'); 
%     scatter(h_ax(ax_idx), x1(~idx_logic),y2(~idx_logic),20,[1 1 1].*0.0,'Marker','s','MarkerFaceColor','w'); 
%     scatter(h_ax(ax_idx), x1(idx_logic),y1(idx_logic),20,[1 1 1].*0.5,"filled"); 
%     scatter(h_ax(ax_idx), x1(idx_logic),y2(idx_logic),20,[1 1 1].*0.0,"filled"); 

    plot(h_ax(ax_idx), x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    plot(h_ax(ax_idx), x2,y4,'-','Color',color_spd{jj},'LineWidth',2);
    
    xlim(h_ax(ax_idx),[0,5]); ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    uistack(h1,'bottom');
%     uistack(h4,'bottom');
    title(h_ax(ax_idx),subject_list(ii)); 
    xlabel(h_ax(ax_idx),'Force (N)');
    
    if ii == 1 && kk == length(nf_target)
        ylabel(h_ax(ax_idx),'Max strain (%)');
    end
end
end

%%
% MAX CURVATURE WITH FORCE
% graph 3 

groupvars = ["Subject","NF_formatted","SPDcond","TrialNbr"]; 
whichstats = {dummyfun1,dummyfun2,dummyfun3,dummyfun4,dummyfun5}; 
targetvar = ["Curv"]; 
tbl = table_normal_loading;
tblstats1 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats1(1:5,:)

groupvars = ["Subject","NF_formatted","SPDcond"]; 
targetvar = ["Fun3_Curv"]; 
tbl = tblstats1;
tblstats2 = grpstats(tbl,groupvars,whichstats,"DataVars",targetvar); 
tblstats2(1:5,:)

ylim_set = [0,1e-2];
title_fig = sprintf("MaxCurv_vs_force"); 
fig = newfig(title_fig); 
set(fig,'Unit','Normalized','Position',[0.0734    0.0287+0.8926/2    0.2641    0.8926/2]); 
Nrow = 1; Ncol = length(subject_list);
h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on'); 
set(h_ax,'FontSize',12,'Clipping','off'); 
for ii = 1:length(subject_list)
subject_ii = subject_list(ii); 
for jj = 1:length(spdcond_list)
    spdcond_jj = spdcond_list(jj);
    
    mask_target = tblstats1.SPDcond == spdcond_jj...
        & tblstats1.Subject == subject_ii; 

    x1 = tblstats1.NF_formatted(mask_target);
    y1 = tblstats1.(targetvar(1))(mask_target);  

    mask_target = tblstats2.SPDcond == spdcond_jj...
        & tblstats2.Subject == subject_ii; 

    x2 = tblstats2.NF_formatted(mask_target);
    y3 = tblstats2.('Fun4_'+targetvar(1))(mask_target); 
    y4 = tblstats2.('Fun5_'+targetvar(1))(mask_target); 
 
    %graph 
    ax_idx = ii; 
    h1 = yline(h_ax(ax_idx),0);
    scatter(h_ax(ax_idx), x1,y1,20,color_spd{jj},'Marker','s'); 
    
    plot(h_ax(ax_idx),x2,y3,'-','Color',color_spd{jj},'LineWidth',2);
    h2 = ciplot(y3-y4,y3+y4,x2,[1 1 1].*0.8,h_ax(ax_idx),0.99);
    
    xlim(h_ax(ax_idx),[0,5]); ylim(h_ax(ax_idx),ylim_set); 
    xlim_get = get(h_ax(ax_idx),'XLim');ylim_get = get(h_ax(ax_idx),'YLim');
    uistack(h1,'bottom');
    title(h_ax(ax_idx),subject_list(ii)); 
    xlabel(h_ax(ax_idx),'Force (N)');
    
    uistack(h3,'bottom');
    uistack(h2,'bottom');
    uistack(h1,'bottom');
    
    if ii == 1 && kk == length(nf_target)
        ylabel(h_ax(ax_idx),'Max curvature (-)');
    end
end
end

