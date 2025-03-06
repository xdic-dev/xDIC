
%Path 
baseDataPath = "E:\transfert_SMP"; 
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 
%Phase 
phase = "slide1";

clear trial subject material 

subject_idx = 0; 
%S01
subject_idx = subject_idx+1; 
subject(subject_idx) = "S01"; 
trial{subject_idx} = [3,4,6,8,9,11,16,20,21,24];%,32,43
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S02
subject_idx = subject_idx+1;
subject(subject_idx) = "S02";
trial{subject_idx} = [1,5,6,8,13,14,16,28,31,48];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S03
subject_idx = subject_idx+1;
subject(subject_idx) = "S03";
trial{subject_idx} = [1,2,3,4,5,7,9,11,13,18,19,26];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S04
subject_idx = subject_idx+1; 
subject(subject_idx) = "S04"; 
trial{subject_idx} = [1,2,3,4,5,6,7,9,12,18,19];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S05
subject_idx = subject_idx+1; 
subject(subject_idx) = "S05"; 
trial{subject_idx} = [2,10,14];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));

% import data
data_multidic = load_DIC_results(...
    'baseDataPath',baseDataPath,...
    'baseResultPath',baseResultPath,...
    'subject',subject,...
    'material',material,...
    'phase',phase,...
    'trial',trial,...
    'deftype','cum'...
    ); 
Nsubject = subject_idx; 
%%
data_multidic([end-1,end],:) = []; 

Ntrial = height(data_multidic);

%% Filtering 
smoothPar.n = 5; 
smoothPar.lambda = 0.5;
smoothspace = 0; 
for itrial = 1:Ntrial
    datai = data_multidic.dic(itrial);
    %J
    [data_multidic.dic(itrial).Jfilt,data_multidic.dic(itrial).Jratefilt] ...
        = filter3Ddeform_time(datai.J,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',datai.Faces,...
        'P',datai.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar); 
    %Epc1
    [data_multidic.dic(itrial).Epc1filt,data_multidic.dic(itrial).Epc1ratefilt] = ...
        filter3Ddeform_time(datai.Epc1,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',datai.Faces,...
        'P',datai.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar);
    %Epc2
    [data_multidic.dic(itrial).Epc2filt,data_multidic.dic(itrial).Epc2ratefilt] = ...
        filter3Ddeform_time(datai.Epc2,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',datai.Faces,...
        'P',datai.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar); 
    %EShearMax
    [data_multidic.dic(itrial).EShearMaxfilt,data_multidic.dic(itrial).EShearMaxratefilt] = ...
        filter3Ddeform_time(datai.EShearMax,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',datai.Faces,...
        'P',datai.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar); 
    fprintf("trial %d done\n",itrial); 
end

%% 
ptile = 5; 
corrcoef_threshold = 0.75; 
ratedeform_threshold = 0.025;
deform_threshold = 0.2; 
threshold = -1; 
results = cell(Ntrial,1); 
for itrial = 1:Ntrial
    %current data 
    subject = data_multidic.subject(itrial);
    trialnumber = data_multidic.trial(itrial);
    dir = data_multidic.dir(itrial);
    nf = data_multidic.nf(itrial);
    spddxl = data_multidic.spddxl(itrial);
    datai=data_multidic.dic(itrial);
    
    Fcentroidlast=datai.FaceCentroids{end};
    Nframe = length(datai.FaceCentroids); 
    
    %refine threshold from last map 
    idx_IN = Fcentroidlast(:,3) < threshold; 
    DM = [Fcentroidlast(idx_IN,1),Fcentroidlast(idx_IN,2),ones(length(Fcentroidlast(idx_IN,1)),1)]; % Design Matrix
    param = DM\Fcentroidlast(idx_IN,3);                   % Estimate Parameters
    B = param; 
    planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.13);
    plane_id = planefit(Fcentroidlast(:,1),Fcentroidlast(:,2)) > Fcentroidlast(:,3); 
    threshold = mean(Fcentroidlast(plane_id,3));
    
    %init
    subject = repelem(subject,Nframe)'; 
    trial = repelem(itrial,Nframe)'; 
    frame = (1:Nframe)'; 
    dir = repelem(dir,Nframe)'; 
    nf = repelem(nf,Nframe)'; 
    spddxl = repelem(spddxl,Nframe)';
    
    z = NaN(Nframe,4); 
    Epc1max = z; Epc2max = z; Jmax = z; Shearmax = z; 
    Epc1min = z; Epc2min = z; Jmin = z; Shearmin = z; 
    Epc1ratemax = z; Epc2ratemax = z; Jratemax = z; Shearratemax = z; 
    Epc1ratemin = z; Epc2ratemin = z; Jratemin = z; Shearratemin = z; 
    Epc1mean = z; Epc2mean = z; Jmean = z; Shearmean = z; 
    
    for ii = 1:Nframe
        Fcentroidnow=datai.FaceCentroids{ii};
        CorrCoef = datai.FaceCorrComb{ii};
        
        [mask_R0,mask_R1,mask_R2,mask_R3] = find_fingerside(...
            'data',datai,...
            'currentFrame',ii);
        
        mask1 = false(length(datai.Epc1filt{ii}),1); 
        mask2 = mask1; 
        
        if length(find(mask_R2==1)) < 200
            fprintf("subject %s - trial %d - dir %s - frame %d suppressed\n",subject{ii},itrial,dir(ii),ii); 
            mask_R2 = false(length(mask_R2),1); 
        end
        if ii~=1
            def = datai.Epc1filt{ii}; 
            defprevious = datai.Epc1filt{ii-1};
            mask1 = abs(def-defprevious)>ratedeform_threshold;
            mask2 = abs(def-1)>deform_threshold;
        end
        corrcoefmask = CorrCoef > corrcoef_threshold;  
        maskDelete = corrcoefmask&mask1&mask2; 
        
        %maximum 
        Epc1max(ii,1)  = prctile(datai.Epc1filt{ii}(mask_R0&~maskDelete),100-ptile); 
        Epc2max(ii,1)  = prctile(datai.Epc2filt{ii}(mask_R0&~maskDelete),100-ptile); 
        Jmax(ii,1)     = prctile(datai.Jfilt{ii}(mask_R0&~maskDelete),100-ptile); 
        Shearmax(ii,1) = prctile(datai.EShearMaxfilt{ii}(mask_R0&~maskDelete),100-ptile); 

        Epc1max(ii,2)  = prctile(datai.Epc1filt{ii}(mask_R1&~maskDelete),100-ptile);
        Epc2max(ii,2)  = prctile(datai.Epc2filt{ii}(mask_R1&~maskDelete),100-ptile);
        Jmax(ii,2)     = prctile(datai.Jfilt{ii}(mask_R1&~maskDelete),100-ptile);
        Shearmax(ii,2) = prctile(datai.EShearMaxfilt{ii}(mask_R1&~maskDelete),100-ptile);

        Epc1max(ii,3)  = prctile(datai.Epc1filt{ii}(mask_R2&~maskDelete),100-ptile);
        Epc2max(ii,3)  = prctile(datai.Epc2filt{ii}(mask_R2&~maskDelete),100-ptile);
        Jmax(ii,3)     = prctile(datai.Jfilt{ii}(mask_R2&~maskDelete),100-ptile);
        Shearmax(ii,3) = prctile(datai.EShearMaxfilt{ii}(mask_R2&~maskDelete),100-ptile);

        Epc1max(ii,4)  = prctile(datai.Epc1filt{ii}(mask_R3&~maskDelete),100-ptile);
        Epc2max(ii,4)  = prctile(datai.Epc2filt{ii}(mask_R3&~maskDelete),100-ptile);
        Jmax(ii,4)     = prctile(datai.Jfilt{ii}(mask_R3&~maskDelete),100-ptile);
        Shearmax(ii,4) = prctile(datai.EShearMaxfilt{ii}(mask_R3&~maskDelete),100-ptile);


        Epc1mean(ii,4)  = mean(datai.Epc1filt{ii}(mask_R3&~maskDelete),'omitnan');
        Epc2mean(ii,4)  = mean(datai.Epc2filt{ii}(mask_R3&~maskDelete),'omitnan');
        Jmean(ii,4)     = mean(datai.Jfilt{ii}(mask_R3&~maskDelete),'omitnan');
        Shearmean(ii,4) = mean(datai.EShearMaxfilt{ii}(mask_R3&~maskDelete),'omitnan');

        %minimum 
        Epc1min(ii,1)  = prctile(datai.Epc1filt{ii}(mask_R0&~maskDelete),ptile); 
        Epc2min(ii,1)  = prctile(datai.Epc2filt{ii}(mask_R0&~maskDelete),ptile); 
        Jmin(ii,1)     = prctile(datai.Jfilt{ii}(mask_R0&~maskDelete),ptile); 
        Shearmin(ii,1) = prctile(datai.EShearMaxfilt{ii}(mask_R0&~maskDelete),ptile); 

        Epc1min(ii,2)  = prctile(datai.Epc1filt{ii}(mask_R1&~maskDelete),ptile);
        Epc2min(ii,2)  = prctile(datai.Epc2filt{ii}(mask_R1&~maskDelete),ptile);
        Jmin(ii,2)     = prctile(datai.Jfilt{ii}(mask_R1&~maskDelete),ptile);
        Shearmin(ii,2) = prctile(datai.EShearMaxfilt{ii}(mask_R1&~maskDelete),ptile);

        Epc1min(ii,3)  = prctile(datai.Epc1filt{ii}(mask_R2&~maskDelete),ptile);
        Epc2min(ii,3)  = prctile(datai.Epc2filt{ii}(mask_R2&~maskDelete),ptile);
        Jmin(ii,3)     = prctile(datai.Jfilt{ii}(mask_R2&~maskDelete),ptile);
        Shearmin(ii,3) = prctile(datai.EShearMaxfilt{ii}(mask_R2&~maskDelete),ptile);

        Epc1min(ii,4)  = prctile(datai.Epc1filt{ii}(mask_R3&~maskDelete),ptile);
        Epc2min(ii,4)  = prctile(datai.Epc2filt{ii}(mask_R3&~maskDelete),ptile);
        Jmin(ii,4)     = prctile(datai.Jfilt{ii}(mask_R3&~maskDelete),ptile);
        Shearmin(ii,4) = prctile(datai.EShearMaxfilt{ii}(mask_R3&~maskDelete),ptile);

        %mean
        Epc1mean(ii,1)  = mean(datai.Epc1filt{ii}(mask_R0&~maskDelete),'omitnan'); 
        Epc2mean(ii,1)  = mean(datai.Epc2filt{ii}(mask_R0&~maskDelete),'omitnan'); 
        Jmean(ii,1)     = mean(datai.Jfilt{ii}(mask_R0&~maskDelete),'omitnan'); 
        Shearmean(ii,1) = mean(datai.EShearMaxfilt{ii}(mask_R0&~maskDelete),'omitnan'); 

        Epc1mean(ii,2)  = mean(datai.Epc1filt{ii}(mask_R1&~maskDelete),'omitnan');
        Epc2mean(ii,2)  = mean(datai.Epc2filt{ii}(mask_R1&~maskDelete),'omitnan');
        Jmean(ii,2)     = mean(datai.Jfilt{ii}(mask_R1&~maskDelete),'omitnan');
        Shearmean(ii,2) = mean(datai.EShearMaxfilt{ii}(mask_R1&~maskDelete),'omitnan');

        Epc1mean(ii,3)  = mean(datai.Epc1filt{ii}(mask_R2&~maskDelete),'omitnan');
        Epc2mean(ii,3)  = mean(datai.Epc2filt{ii}(mask_R2&~maskDelete),'omitnan');
        Jmean(ii,3)     = mean(datai.Jfilt{ii}(mask_R2&~maskDelete),'omitnan');
        Shearmean(ii,3) = mean(datai.EShearMaxfilt{ii}(mask_R2&~maskDelete),'omitnan');

        Epc1mean(ii,4)  = prctile(datai.Epc1filt{ii}(mask_R3&~maskDelete),ptile);
        Epc2mean(ii,4)  = prctile(datai.Epc2filt{ii}(mask_R3&~maskDelete),ptile);
        Jmean(ii,4)     = prctile(datai.Jfilt{ii}(mask_R3&~maskDelete),ptile);
        Shearmean(ii,4) = prctile(datai.EShearMaxfilt{ii}(mask_R3&~maskDelete),ptile);


        %rate
        %maximum 
        Epc1ratemax(ii,1)  = prctile(datai.Epc1ratefilt{ii}(mask_R0&~maskDelete),100-ptile); 
        Epc2ratemax(ii,1)  = prctile(datai.Epc2ratefilt{ii}(mask_R0&~maskDelete),100-ptile); 
        Jratemax(ii,1)     = prctile(datai.Jratefilt{ii}(mask_R0&~maskDelete),100-ptile); 
        Shearratemax(ii,1) = prctile(datai.EShearMaxratefilt{ii}(mask_R0&~maskDelete),100-ptile); 

        Epc1ratemax(ii,2)  = prctile(datai.Epc1ratefilt{ii}(mask_R1&~maskDelete),100-ptile);
        Epc2ratemax(ii,2)  = prctile(datai.Epc2ratefilt{ii}(mask_R1&~maskDelete),100-ptile);
        Jratemax(ii,2)     = prctile(datai.Jratefilt{ii}(mask_R1&~maskDelete),100-ptile);
        Shearratemax(ii,2) = prctile(datai.EShearMaxratefilt{ii}(mask_R1&~maskDelete),100-ptile);

        Epc1ratemax(ii,3)  = prctile(datai.Epc1ratefilt{ii}(mask_R2&~maskDelete),100-ptile); 
        Epc2ratemax(ii,3)  = prctile(datai.Epc2ratefilt{ii}(mask_R2&~maskDelete),100-ptile); 
        Jratemax(ii,3)     = prctile(datai.Jratefilt{ii}(mask_R2&~maskDelete),100-ptile); 
        Shearratemax(ii,3) = prctile(datai.EShearMaxratefilt{ii}(mask_R2&~maskDelete),100-ptile); 

        Epc1ratemax(ii,4)  = prctile(datai.Epc1ratefilt{ii}(mask_R3&~maskDelete),100-ptile);
        Epc2ratemax(ii,4)  = prctile(datai.Epc2ratefilt{ii}(mask_R3&~maskDelete),100-ptile);
        Jratemax(ii,4)     = prctile(datai.Jratefilt{ii}(mask_R3&~maskDelete),100-ptile);
        Shearratemax(ii,4) = prctile(datai.EShearMaxratefilt{ii}(mask_R3&~maskDelete),100-ptile);


        %minimum 
        Epc1ratemin(ii,1)  = prctile(datai.Epc1ratefilt{ii}(mask_R0&~maskDelete),ptile); 
        Epc2ratemin(ii,1)  = prctile(datai.Epc2ratefilt{ii}(mask_R0&~maskDelete),ptile); 
        Jratemin(ii,1)     = prctile(datai.Jratefilt{ii}(mask_R0&~maskDelete),ptile); 
        Shearratemin(ii,1) = prctile(datai.EShearMaxratefilt{ii}(mask_R0&~maskDelete),ptile); 

        Epc1ratemin(ii,2)  = prctile(datai.Epc1ratefilt{ii}(mask_R1&~maskDelete),ptile);
        Epc2ratemin(ii,2)  = prctile(datai.Epc2ratefilt{ii}(mask_R1&~maskDelete),ptile);
        Jratemin(ii,2)     = prctile(datai.Jratefilt{ii}(mask_R1&~maskDelete),ptile);
        Shearratemin(ii,2) = prctile(datai.EShearMaxratefilt{ii}(mask_R1&~maskDelete),ptile);

        Epc1ratemin(ii,3)  = prctile(datai.Epc1ratefilt{ii}(mask_R2&~maskDelete),ptile); 
        Epc2ratemin(ii,3)  = prctile(datai.Epc2ratefilt{ii}(mask_R2&~maskDelete),ptile); 
        Jratemin(ii,3)     = prctile(datai.Jratefilt{ii}(mask_R2&~maskDelete),ptile); 
        Shearratemin(ii,3) = prctile(datai.EShearMaxratefilt{ii}(mask_R2&~maskDelete),ptile); 

        Epc1ratemin(ii,4)  = prctile(datai.Epc1ratefilt{ii}(mask_R3&~maskDelete),ptile);
        Epc2ratemin(ii,4)  = prctile(datai.Epc2ratefilt{ii}(mask_R3&~maskDelete),ptile);
        Jratemin(ii,4)     = prctile(datai.Jratefilt{ii}(mask_R3&~maskDelete),ptile);
        Shearratemin(ii,4) = prctile(datai.EShearMaxratefilt{ii}(mask_R3&~maskDelete),ptile);
        
    end
    fprintf("trial %d done\n",itrial); 
    results{itrial} = table(subject,trial,frame,dir,nf,spddxl,...
        Epc1max,Epc2max,Jmax,Shearmax,...
        Epc1min,Epc2min,Jmin,Shearmin,...
        Epc1ratemax,Epc2ratemax,Jratemax,Shearratemax,...
        Epc1ratemin,Epc2ratemin,Jratemin,Shearratemin,...
        Epc1mean,Epc2mean,Jmean,Shearmean); 
%     results = [results;result];
end
%%

colorStretch1 = rgb('LightSalmon');
colorStretch2 = rgb('FireBrick');
colorCompr1 = rgb('DeepSkyBlue');
colorCompr2 = rgb('DarkBlue');

colorR{1} = rgb('BurlyWood'); 
colorR{2} = rgb('LightCoral'); 
colorR{3} = rgb('CadetBlue'); 
colorR{4} = rgb('PaleGreen'); 

ratefactor = 25; 
dirvec = ["Ubnf","Rbnf","Pbnf","Dbnf"]; 
deform1 = "Epc1"; 
deform2 = "Epc2"; 
subjectlist = unique(data_multidic.subject); 
%results_cum = results; 
title_fig = "first and second principal strain all subjects"; 
fig = newfig(title_fig); 
ax = subplot_ax(4,4); hold(ax,'on'); 
set(ax,'Color','None'); 
msize = 0; 

plot(ax(4),[-45*3.8,45],[-20,-20;0,0;20,20],'-','Color',[1 1 1].*0.7)
plot(ax(8),[-45*3.8,45],[-20,-20;0,0;20,20],'-','Color',[1 1 1].*0.7)
plot(ax(12),[-45*3.8,45],[-20,-20;0,0;20,20],'-','Color',[1 1 1].*0.7)
plot(ax(16),[-45*3.8,45],[-20,-20;0,0;20,20],'-','Color',[1 1 1].*0.7)

% subjectNum = find(strcmp(data_multidic_cum.subject(itrial),subjectlist));
for iregion = 1:4
    for itrial = 1:Ntrial
        ii = 0; 
        for idir = dirvec
            ii = ii+1; 
            ax_idx = ii+(iregion-1)*4;
            colorBG = colorR{ceil(ax_idx/4)}; 
            if data_multidic.spd(itrial)==5&&strcmp(data_multidic.dir(itrial),idir)%&&data_multidic.nf(itrial)==5
                if data_multidic.nf(itrial)==1
                    coloritrial = [1 1 1].*0.5;
                    plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results{itrial}.(deform1+"ratemin")(:,iregion)*1e2*ratefactor,'.-','Color',colorCompr2,'MarkerSize',msize+5)
                    plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results{itrial}.(deform2+"ratemax")(:,iregion)*1e2*ratefactor,'.-','Color',colorStretch2,'MarkerSize',msize+5)
                else 
                    coloritrial = [1 1 1].*0.0;
                    plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results{itrial}.(deform1+"ratemin")(:,iregion)*1e2*ratefactor,'.-','Color',colorCompr2,'MarkerSize',msize+5)
                    plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results{itrial}.(deform2+"ratemax")(:,iregion)*1e2*ratefactor,'.-','Color',colorStretch2,'MarkerSize',msize+5)
                end
                ylim(ax(ax_idx),[-30,30]); 
            else
            end
            %set(ax(ax_idx),'Color',colorBG); 
            ax(ax_idx).XAxis.Visible = 'off';
            ax(ax_idx).YAxis.Visible = 'off';
            ax(ax_idx).XTick = [];
            ax(ax_idx).YTick = [];
        end
    end
end
set(ax,'Clipping','off');
set(ax,'FontSize',15); 
linkaxes(ax,'x')
xlim(ax(end),[1 45]); 

title(ax(1),"Ulnar"); 
title(ax(2),"Radial"); 
title(ax(3),"Proximal");
title(ax(4),"Distal");

% ax(end).XAxis.Visible = 'on';
% ax(end).XTick = [1,20,40];
% ax(end).XTickLabel = num2str(ax(ax_idx).XTick');
% ax(end-1).XAxis.Visible = 'on';
% ax(end-1).XTick = [1,20,40];
% ax(end-1).XTickLabel = num2str(ax(ax_idx).XTick');
% ax(end-2).XAxis.Visible = 'on'; 
% ax(end-2).XTick = [1,20,40];
% ax(end-2).XTickLabel = num2str(ax(ax_idx).XTick');
% ax(end-3).XAxis.Visible = 'on';
% ax(end-3).XTick = [1,20,40];
% ax(end-3).XTickLabel = num2str(ax(ax_idx).XTick');

ylabel(ax(1),"R0",'Color',colorR{1},'FontWeight','bold'); 
ax(1).YAxis.Visible = 'on';
ax(1).YTick = [-40,-20,0,20,40];
ax(1).YTickLabel = num2str(ax(1).YTick');
ylabel(ax(5),"R1",'Color',colorR{2},'FontWeight','bold'); 
ax(5).YAxis.Visible = 'on';
ax(5).YTick = [-40,-20,0,20,40];
ax(5).YTickLabel = num2str(ax(5).YTick');
ylabel(ax(9),"R2",'Color',colorR{3},'FontWeight','bold'); 
ax(9).YAxis.Visible = 'on';
ax(9).YTick = [-40,-20,0,20,40];
ax(9).YTickLabel = num2str(ax(9).YTick');
ylabel(ax(13),"R3",'Color',colorR{4},'FontWeight','bold'); 
ax(13).YAxis.Visible = 'on';
ax(13).YTick = [-40,-20,0,20,40];
ax(13).YTickLabel = num2str(ax(13).YTick');

%total 
text(ax(1),0,40,["Total (%)"],...%["Total","%"]
        'FontSize',12,...
        'Rotation',0,...
        'FontWeight','bold',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle');

%Xaxis 
y0 = -30;
plot(ax(end),[45 (45-25)],[y0 y0],'k-','linew',2);
text(ax(end),mean([45 (45-25)]),y0-1,"1000ms",...
        'FontSize',15,...
        'Rotation',0,...
        'FontWeight','bold',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top');
plot(ax(12),[45 (45-25)],[y0 y0],'k-','linew',2);
text(ax(12),mean([45 (45-25)]),y0-1,"1000ms",...
    'FontSize',15,...
    'Rotation',0,...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
plot(ax(8),[45 (45-25)],[y0 y0],'k-','linew',2);
text(ax(8),mean([45 (45-25)]),y0-1,"1000ms",...
    'FontSize',15,...
    'Rotation',0,...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
plot(ax(4),[45 (45-25)],[y0 y0],'k-','linew',2);
text(ax(4),mean([45 (45-25)]),y0-1,"1000ms",...
    'FontSize',15,...
    'Rotation',0,...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
    
%legend
qw = cell(2,1); 
qw{1} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],colorCompr2);
qw{2} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],colorStretch2);
% qw{3} = plot(ax(ax_idx),NaN,NaN,'-','Color','k');
% qw{4} = plot(ax(ax_idx),NaN,NaN,'--','Color','k');
if strcmp(deform1,"J")
legendlabel = {'Min surf. ch.','Max surf. ch.'}; 
else
legendlabel = {'Max COMPR.','Max STRET.'}; 
end  
hlegend = legend(ax(7),[qw{:}],legendlabel);%'min','max' 
set(hlegend, 'NumColumns', 2,...
    'Location','best',...
    'Box','off',...
    'FontSize',15);
box(hlegend,'off'); 




















