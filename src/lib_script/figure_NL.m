
%Path 
baseDataPath = "E:\transfert_SMP"; 
baseResultPath = "D:\INMACOSY Dropbox\Donatien Doumont\myDropbox\SMP_data\MultiView_project\2023_04_acquisition_sept_oct"; 
%Phase 
phase = "loading";
%Material
% material = "coating"; 

subject_idx = 0; 
%S01
subject_idx = subject_idx+1; 
subject(subject_idx) = "S01"; 
trial{subject_idx} = [9,11,45,46];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S02
subject_idx = subject_idx+1;
subject(subject_idx) = "S02";
trial{subject_idx} = [1,2,47,68];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S03
subject_idx = subject_idx+1;
subject(subject_idx) = "S03";
trial{subject_idx} = [1,6];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S04
subject_idx = subject_idx+1; 
subject(subject_idx) = "S04"; 
trial{subject_idx} = [1,6,25,26];
material{subject_idx} = repelem("coating",length(trial{subject_idx}));
%S05
subject_idx = subject_idx+1; 
subject(subject_idx) = "S05"; 
trial{subject_idx} = [10,15,28,48];
material{subject_idx} = repelem("glass",length(trial{subject_idx}));
%S06
subject_idx = subject_idx+1; 
subject(subject_idx) = "S06"; 
trial{subject_idx} = [7,11,26,30];
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
data_multidic([end-1,end],:) = []; 
Ntrial = height(data_multidic);

%% Filtering 
smoothPar.n = 5; 
smoothPar.lambda = 0.5;
smoothspace = 0; 
for itrial = 1:Ntrial
    data = data_multidic.dic(itrial);
    %J
    [data_multidic.dic(itrial).Jfilt,data_multidic.dic(itrial).Jratefilt] ...
        = filter3Ddeform_time(data.J-1,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',data.Faces,...
        'P',data.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar); 
    %Epc1
    [data_multidic.dic(itrial).Epc1filt,data_multidic.dic(itrial).Epc1ratefilt] = ...
        filter3Ddeform_time(data.Epc1,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',data.Faces,...
        'P',data.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar);
    %Epc2
    [data_multidic.dic(itrial).Epc2filt,data_multidic.dic(itrial).Epc2ratefilt] = ...
        filter3Ddeform_time(data.Epc2,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',data.Faces,...
        'P',data.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar); 
    %EShearMax
    [data_multidic.dic(itrial).EShearMaxfilt,data_multidic.dic(itrial).EShearMaxratefilt] = ...
        filter3Ddeform_time(data.EShearMax,...
        'freqFilt',0.07,...
        'freqAcq',1,...
        'F',data.Faces,...
        'P',data.Points3D,...
        'smoothSpace',smoothspace,...
        'smoothPar',smoothPar); 
end
%%
% title_fig = 'test'; fig = newfig(title_fig); ax = subplot_ax(2,1);
% hold(ax,'on'); fps = 50; ax_idx = 1; ptile = 5; mask =
% strcmp(data_multidic.dir,"Ubnf"); mask = mask & data_multidic.nf==1; for
% itrial = find(mask==1)'
%     data = data_multidic.dic(itrial).Deform; Nframe = length(data); val1
%     = zeros(1,Nframe); val2 = zeros(1,Nframe); for ii = 1:Nframe
%         val1(ii) = prctile(data{ii},100-ptile)-1; val2(ii) =
%         -(prctile(data{ii},ptile)-1);
%     end max(val1) max(val2) if data_multidic.spd(itrial)==5
%         displ =
%         linspace(0,data_multidic.spd(itrial)*Nframe/fps*2,Nframe); va11 =
%         (val1)*100*fps/2; val2 = (val2)*100*fps/2;
%     else
%         displ =
%         linspace(0,data_multidic.spd(itrial)*Nframe/fps*1,Nframe); va11 =
%         (val1)*100*fps; val2 = (val2)*100*fps;
%     end if strcmp(data_multidic.subject(itrial),'S01')
%         plot(ax(ax_idx),displ,val1,'r--');
%         plot(ax(ax_idx),displ,val2,'b--');
%     elseif strcmp(data_multidic.subject(itrial),'S03')
%         plot(ax(ax_idx),displ,val1,'r-.');
%         plot(ax(ax_idx),displ,val2,'b-.');
%     else
%         plot(ax(ax_idx),displ,val1,'r-');
%         plot(ax(ax_idx),displ,val2,'b-');
%     end
% end ax(ax_idx).XLabel.String = "displ [mm]"; ax(ax_idx).YLabel.String =
% "strain rate [%/s]";
% 
% ax_idx = 2; ptile = 5; mask = strcmp(data_multidic.dir,"Ubnf"); mask =
% mask & data_multidic.nf==5; for itrial = find(mask==1)'
%     data = data_multidic.dic(itrial).Deform; Nframe = length(data); val1
%     = zeros(1,Nframe); val2 = zeros(1,Nframe); if
%     data_multidic.spd(itrial)==5
%         displ = linspace(0,Nframe/fps*2,Nframe);
%     elseif data_multidic.spd(itrial)==10
%         displ = linspace(0,Nframe/fps*1,Nframe);
%     end
%     
%     for ii = 1:Nframe
%         val1(ii) = prctile(data{ii},100-ptile)-1; val2(ii) =
%         -(prctile(data{ii},ptile)-1);
%     end if strcmp(data_multidic.subject(itrial),'S01')
%         plot(ax(ax_idx),displ,val1,'r--');
%         plot(ax(ax_idx),displ,val2,'b--');
%     elseif strcmp(data_multidic.subject(itrial),'S03')
%         plot(ax(ax_idx),displ,val1,'r-.');
%         plot(ax(ax_idx),displ,val2,'b-.');
%     else
%         plot(ax(ax_idx),displ,val1,'r-');
%         plot(ax(ax_idx),displ,val2,'b-');
%     end
% end
%% Segmentation of data 
a = {NaN,NaN}';
data_multidic = addvars(data_multidic,[repelem(a,1,Ntrial)]','NewVariableNames','results');
%%
for ii = 1:Nsubject
    Filterfrequency = 10;
    datafilepath = fullfile(baseDataPath,"rawdata",subject(ii),'speckles',material{ii}(1),"robot");
    data_robot{ii} = import_robotdata_subject(datafilepath,...
        'FreqFilter',Filterfrequency); 
    data_nf_idx_ii = load(fullfile(baseResultPath,...
        'analysis',subject(ii),'frameidx.mat')); 
    data_nf_idx{ii} = data_nf_idx_ii.frameidx{strcmp(material{ii}(1),"coating")+1}.loading; 
   
end
%%
kk = 0;
for ii = 1:Nsubject
    for itrial = trial{ii}
        kk = kk+1; 
        
        if strcmp(data_multidic.spddxl(kk),"lowspd")
            iframe = data_nf_idx{ii}.start(itrial):2:data_nf_idx{ii}.end(itrial);
            iframe = round(6.5e3-0.5e3+iframe*1e3/25);
        else
            iframe = data_nf_idx{ii}.start(itrial):1:data_nf_idx{ii}.end(itrial);
            iframe = round(6.5e3-0.5e3+iframe*1e3/50);
        end
        nf_frame{kk} = data_robot{ii}{itrial}.nf(iframe); 
    end
end
%%
figure; 
plot(data_robot{end}{itrial}.time, data_robot{end}{itrial}.nf);hold on; 
iframe(end)
plot(data_robot{end}{itrial(end)}.time(iframe(1)), data_robot{end}{itrial(end)}.nf(iframe(1)),'ro');
plot(data_robot{end}{itrial(end)}.time(iframe(end)), data_robot{end}{itrial(end)}.nf(iframe(end)),'ro');
%%
ptile = 5; 
corrcoef_threshold = 0.75; 
ratedeform_threshold = 0.025;
deform_threshold = 0.2; 
threshold = -1; 
results = cell(Ntrial,1); 
for itrial = 1:Ntrial
    subject = data_multidic.subject(itrial);
    trialnumber = data_multidic.trial(itrial);
    nf = data_multidic.nf(itrial);
    spddxl = data_multidic.spddxl(itrial);
    datai=data_multidic.dic(itrial);
    
%     varNames = ["subject","trial","frame","nf","spddxl","Epc1maxIN","Epc2maxIN","JmaxIN","ShearmaxIN",...
%             "Epc1maxOUT","Epc2maxOUT","JmaxOUT","ShearmaxOUT",...
%             "Epc1minIN","Epc2minIN","JminIN","ShearminIN",...
%             "Epc1minOUT","Epc2minOUT","JminOUT","ShearminOUT"];
%     varTypes = ["string","string","double","double","string",repelem("double",16)];
%     result = table('Size',[data_multidic.Nframe(itrial) length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
 
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
    frame = [1:Nframe]'; 
    nf = repelem(nf,Nframe)'; 
    spddxl = repelem(spddxl,Nframe)';
    
    z = zeros(Nframe,2); 
    Epc1max = z; Epc2max = z; Jmax = z; Shearmax = z; 
    Epc1min = z; Epc2min = z; Jmin = z; Shearmin = z; 
    Epc1ratemax = z; Epc2ratemax = z; Jratemax = z; Shearratemax = z; 
    Epc1ratemin = z; Epc2ratemin = z; Jratemin = z; Shearratemin = z; 
    
    for ii = 1:Nframe
        Fcentroidnow=datai.FaceCentroids{ii};
        CorrCoef = datai.FaceCorrComb{ii};
        %fit plane with point below threshold
        IN_threshold = 0; %200
        idx_IN = Fcentroidnow(:,3) < threshold+0.15; 
        DM = [Fcentroidnow(idx_IN,1),Fcentroidnow(idx_IN,2),ones(length(Fcentroidnow(idx_IN,1)),1)]; % Design Matrix
        param = DM\Fcentroidnow(idx_IN,3);                   % Estimate Parameters
        B = param; 
        planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.15);

        %find the boundary 
        plane_id = Fcentroidnow(:,3) < planefit(Fcentroidnow(:,1),Fcentroidnow(:,2)); 
        if length(find(idx_IN==1))<IN_threshold
            plane_id = false(length(datai.Epc1filt{ii}),1); 
        end
        %
        mask1 = false(length(datai.Epc1filt{ii}),1); 
        mask2 = mask1; 
        
        if ii~=1
            def = datai.Epc1filt{ii}; 
            defprevious = datai.Epc1filt{ii-1};
            mask1 = abs(def-defprevious)>ratedeform_threshold;
            mask2 = abs(def-1)>deform_threshold;
        end
        corrcoefmask = CorrCoef > corrcoef_threshold;  
        maskDelete = corrcoefmask&mask1&mask2; 
        
        if true
            %cum
            %maximum 
            Epc1max(ii,1)  = prctile(datai.Epc1filt{ii}(plane_id&~maskDelete),100-ptile); 
            Epc2max(ii,1)  = prctile(datai.Epc2filt{ii}(plane_id&~maskDelete),100-ptile); 
            Jmax(ii,1)     = prctile(datai.Jfilt{ii}(plane_id&~maskDelete),100-ptile); 
            Shearmax(ii,1) = prctile(datai.EShearMaxfilt{ii}(plane_id&~maskDelete),100-ptile); 

            Epc1max(ii,2)  = prctile(datai.Epc1filt{ii}(~plane_id&~maskDelete),100-ptile);
            Epc2max(ii,2)  = prctile(datai.Epc2filt{ii}(~plane_id&~maskDelete),100-ptile);
            Jmax(ii,2)     = prctile(datai.Jfilt{ii}(~plane_id&~maskDelete),100-ptile);
            Shearmax(ii,2) = prctile(datai.EShearMaxfilt{ii}(~plane_id&~maskDelete),100-ptile);

            %minimum 
            Epc1min(ii,1)  = prctile(datai.Epc1filt{ii}(plane_id&~maskDelete),ptile); 
            Epc2min(ii,1)  = prctile(datai.Epc2filt{ii}(plane_id&~maskDelete),ptile); 
            Jmin(ii,1)     = prctile(datai.Jfilt{ii}(plane_id&~maskDelete),ptile); 
            Shearmin(ii,1) = prctile(datai.EShearMaxfilt{ii}(plane_id&~maskDelete),ptile); 

            Epc1min(ii,2)  = prctile(datai.Epc1filt{ii}(~plane_id&~maskDelete),ptile);
            Epc2min(ii,2)  = prctile(datai.Epc2filt{ii}(~plane_id&~maskDelete),ptile);
            Jmin(ii,2)     = prctile(datai.Jfilt{ii}(~plane_id&~maskDelete),ptile);
            Shearmin(ii,2) = prctile(datai.EShearMaxfilt{ii}(~plane_id&~maskDelete),ptile);
            
            %rate
            %maximum 
            Epc1ratemax(ii,1)  = prctile(datai.Epc1ratefilt{ii}(plane_id&~maskDelete),100-ptile); 
            Epc2ratemax(ii,1)  = prctile(datai.Epc2ratefilt{ii}(plane_id&~maskDelete),100-ptile); 
            Jratemax(ii,1)     = prctile(datai.Jratefilt{ii}(plane_id&~maskDelete),100-ptile); 
            Shearratemax(ii,1) = prctile(datai.EShearMaxratefilt{ii}(plane_id&~maskDelete),100-ptile); 

            Epc1ratemax(ii,2)  = prctile(datai.Epc1ratefilt{ii}(~plane_id&~maskDelete),100-ptile);
            Epc2ratemax(ii,2)  = prctile(datai.Epc2ratefilt{ii}(~plane_id&~maskDelete),100-ptile);
            Jratemax(ii,2)     = prctile(datai.Jratefilt{ii}(~plane_id&~maskDelete),100-ptile);
            Shearratemax(ii,2) = prctile(datai.EShearMaxratefilt{ii}(~plane_id&~maskDelete),100-ptile);

            %minimum 
            Epc1ratemin(ii,1)  = prctile(datai.Epc1ratefilt{ii}(plane_id&~maskDelete),ptile); 
            Epc2ratemin(ii,1)  = prctile(datai.Epc2ratefilt{ii}(plane_id&~maskDelete),ptile); 
            Jratemin(ii,1)     = prctile(datai.Jratefilt{ii}(plane_id&~maskDelete),ptile); 
            Shearratemin(ii,1) = prctile(datai.EShearMaxratefilt{ii}(plane_id&~maskDelete),ptile); 

            Epc1ratemin(ii,2)  = prctile(datai.Epc1ratefilt{ii}(~plane_id&~maskDelete),ptile);
            Epc2ratemin(ii,2)  = prctile(datai.Epc2ratefilt{ii}(~plane_id&~maskDelete),ptile);
            Jratemin(ii,2)     = prctile(datai.Jratefilt{ii}(~plane_id&~maskDelete),ptile);
            Shearratemin(ii,2) = prctile(datai.EShearMaxratefilt{ii}(~plane_id&~maskDelete),ptile);
        else 
            %maximum 
            Epc1max(ii,1)  = prctile(datai.Epc1{ii}(plane_id&~maskDelete),100-ptile); 
            Epc2max(ii,1)  = prctile(datai.Epc2{ii}(plane_id&~maskDelete),100-ptile); 
            Jmax(ii,1)     = prctile(datai.J{ii}(plane_id&~maskDelete),100-ptile); 
            Shearmax(ii,1) = prctile(datai.EShearMax{ii}(plane_id&~maskDelete),100-ptile); 

            Epc1max(ii,2)  = prctile(datai.Epc1{ii}(~plane_id&~maskDelete),100-ptile);
            Epc2max(ii,2)  = prctile(datai.Epc2{ii}(~plane_id&~maskDelete),100-ptile);
            Jmax(ii,2)     = prctile(datai.J{ii}(~plane_id&~maskDelete),100-ptile);
            Shearmax(ii,2) = prctile(datai.EShearMax{ii}(~plane_id&~maskDelete),100-ptile);

            %minimum 
            Epc1min(ii,1)  = prctile(datai.Epc1{ii}(plane_id&~maskDelete),ptile); 
            Epc2min(ii,1)  = prctile(datai.Epc2{ii}(plane_id&~maskDelete),ptile); 
            Jmin(ii,1)     = prctile(datai.J{ii}(plane_id&~maskDelete),ptile); 
            Shearmin(ii,1) = prctile(datai.EShearMax{ii}(plane_id&~maskDelete),ptile); 

            Epc1min(ii,2)  = prctile(datai.Epc1{ii}(~plane_id&~maskDelete),ptile);
            Epc2min(ii,2)  = prctile(datai.Epc2{ii}(~plane_id&~maskDelete),ptile);
            Jmin(ii,2)     = prctile(datai.J{ii}(~plane_id&~maskDelete),ptile);
            Shearmin(ii,2) = prctile(datai.EShearMax{ii}(~plane_id&~maskDelete),ptile);
        end
    end
    results{itrial} = table(subject,trial,frame,nf,spddxl,...
        Epc1max,Epc2max,Jmax,Shearmax,...
        Epc1min,Epc2min,Jmin,Shearmin,...
        Epc1ratemax,Epc2ratemax,Jratemax,Shearratemax,...
        Epc1ratemin,Epc2ratemin,Jratemin,Shearratemin); 
%     results = [results;result];
end
%% for individual subject 
colorIN = rgb('Chocolate');
colorOUT = rgb('DarkSlateGray');

deform1 = "Epc1"; 
deform2 = "Epc2"; 
subjectlist = unique(data_multidic.subject); 
%results_cum = results; 
title_fig = "cum : first and second principal strain all subjects"; 
fig = newfig(title_fig); 
ax = subplot_ax(length(subjectlist),4); hold(ax,'on'); 
hold on; 
data_multidic.Nframe(1)
size(results_cum{1}.(deform1+"min")(:,1))
for itrial = 1:Ntrial
    subjectNum = find(strcmp(data_multidic_cum.subject(itrial),subjectlist));
    if data_multidic.nf(itrial)==1&&strcmp(data_multidic.spddxl(itrial),"lowspd")
        ax_idx = 1+(subjectNum-1)*4; 
        length(nf_frame{itrial})
        length(results_cum{itrial}.(deform1+"min")(:,1))
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JmaxIN,'y')
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JminOUT,'g')
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"minOUT"),'--','Color',colorOUT)
        plot(ax(ax_idx),nf_frame{itrial},results_cum{itrial}.(deform1+"min")(:,1),'--','Color',colorIN)
        plot(ax(ax_idx),nf_frame{itrial},results_cum{itrial}.(deform2+"max")(:,2),'-','Color',colorOUT)
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"maxIN"),'-','Color',colorIN)
        xlim(ax(ax_idx),[0,5]); 
        ylim(ax(ax_idx),[-0.15,0.15]);
        title(ax(ax_idx),subjectlist(subjectNum)+"-1N-lowspd"); 
    elseif data_multidic.nf(itrial)==5&&strcmp(data_multidic.spddxl(itrial),"lowspd")
        ax_idx = 3+(subjectNum-1)*4; 
        
        length(nf_frame{itrial})
        length(results_cum{itrial}.(deform1+"min")(:,1))
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JmaxIN,'y')
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JminOUT,'g')
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"minOUT"),'--','Color',colorOUT)
        plot(ax(ax_idx),nf_frame{itrial},results_cum{itrial}.(deform1+"min")(:,1),'--','Color',colorIN)
        plot(ax(ax_idx),nf_frame{itrial},results_cum{itrial}.(deform2+"max")(:,2),'-','Color',colorOUT)
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"maxIN"),'-','Color',colorIN)
        xlim(ax(ax_idx),[0,5]); 
        ylim(ax(ax_idx),[-0.15,0.15]);
        title(ax(ax_idx),subjectlist(subjectNum)+"-5N-lowspd");  
    end
    
    if data_multidic.nf(itrial)==1&&strcmp(data_multidic.spddxl(itrial),"highspd")
        ax_idx = 2+(subjectNum-1)*4; 
        length(nf_frame{itrial})
        length(results_cum{itrial}.(deform1+"min")(:,1))
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JmaxIN,'y')
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JminOUT,'g')
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"minOUT"),'--','Color',colorOUT)
        plot(ax(ax_idx),nf_frame{itrial},results_cum{itrial}.(deform1+"min")(:,1),'--','Color',colorIN)
        plot(ax(ax_idx),nf_frame{itrial},results_cum{itrial}.(deform2+"max")(:,2),'-','Color',colorOUT)
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"maxIN"),'-','Color',colorIN)
        xlim(ax(ax_idx),[0,5]); 
        ylim(ax(ax_idx),[-0.15,0.15]); 
        title(ax(ax_idx),subjectlist(subjectNum)+"-1N-highspd"); 
    elseif data_multidic.nf(itrial)==5&&strcmp(data_multidic.spddxl(itrial),"highspd")
        ax_idx = 4+(subjectNum-1)*4; 
        length(nf_frame{itrial})
        length(results_cum{itrial}.(deform1+"min")(:,1))
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JmaxIN,'y')
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JminOUT,'g')
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"minOUT"),'--','Color',colorOUT)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"min")(:,1),'--','Color',colorIN)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"max")(:,2),'-','Color',colorOUT)
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"maxIN"),'-','Color',colorIN)
        xlim(ax(ax_idx),[0,5]); 
        ylim(ax(ax_idx),[-0.15,0.15]); 
        title(ax(ax_idx),subjectlist(subjectNum)+"-5N-highspd"); 
    end
end
ax_idx = 4; 
title(ax(ax_idx),'5N-highspd'); 
%legend
qw = cell(4,1); 
qw{1} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],colorOUT);
qw{2} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],colorIN);
qw{3} = plot(ax(ax_idx),NaN,NaN,'-','Color','k');
qw{4} = plot(ax(ax_idx),NaN,NaN,'--','Color','k');
hlegend = legend(ax(ax_idx),[qw{:}],{'IN','OUT','min','max'}); 
set(hlegend, 'NumColumns', 2,...
    'Location','best',...
    'Box','off',...
    'FontSize',10);
box(hlegend,'off'); 
%% for all subject 
colorIN = rgb('DarkBlue');
colorOUT = rgb('FireBrick');

ratefactor = 25; 

deform1 = "Epc1"; 
deform2 = "Epc2"; 
subjectlist = unique(data_multidic.subject); 
% results_cum = results; 

title_fig = "cum : first and second principal strain all subjects 075"; 
fig = newfig(title_fig); 
set(fig,'Position',[12,1,18,15]); 
ax = subplot_ax(2,4,'merge',{[2,3],[6,7]},...
    'shownum',0); 
hold(ax,'on'); 
 
msize = 10; 

plot(ax(4),[-10.5,1.5+1],[-10,-10;-5,-5;0,0;5,5;10,10],'-','Color',[1 1 1].*0.7)
plot(ax(8),[-10.5,1.5+1],[-60,-60;-40,-40;-20,-20;0,0;20,20],'-','Color',[1 1 1].*0.7)

mask = data_multidic.nf==1&strcmp(data_multidic.spddxl,"lowspd"); 
cumulativeIN1 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
cumulativeOUT1 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
rateIN1 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
rateOUT1 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
ii = 0; 
for itrial = 1:Ntrial
    subjectNum = find(strcmp(data_multidic.subject(itrial),subjectlist));
    if data_multidic.nf(itrial)==1&&strcmp(data_multidic.spddxl(itrial),"lowspd")
        ii = ii+1; 
        ax_idx = 1;
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform1+"max")(1:end-2,1)*1e2,'y')
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform2+"min")(1:end-2,1)*1e2,'g')
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"minOUT"),'--','Color',colorOUT)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform1+"min")(1:end-2,1)*1e2,'.-','Color',[1 1 1].*0.5,'MarkerSize',msize)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform2+"max")(1:end-2,2)*1e2,'.-','Color',[1 1 1].*0.5,'MarkerSize',msize)
        cumulativeIN1(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform1+"min")(1:end-2,1)*1e2;
        cumulativeOUT1(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform2+"max")(1:end-2,2)*1e2;

%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"maxIN"),'-','Color',colorIN)
        xlim(ax(ax_idx),[1,data_multidic.Nframe(itrial)]); 
        ylim(ax(ax_idx),[-0.1*1e2,0.1*1e2]);

        ax_idx = 5;
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform1+"ratemin")(1:end-2,1)*ratefactor*1e2,'.--','Color',[1 1 1].*0.5,'MarkerSize',msize)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform2+"ratemax")(1:end-2,2)*ratefactor*1e2,'.--','Color',[1 1 1].*0.5,'MarkerSize',msize)
        
        rateIN1(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform1+"ratemin")(1:end-2,1)*ratefactor*1e2;
        rateOUT1(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform2+"ratemax")(1:end-2,2)*ratefactor*1e2;

        xlim(ax(ax_idx),[1,15]); 
        ylim(ax(ax_idx),[-0.03*ratefactor*1e2,0.01*ratefactor*1e2]);

    end
end
%title 
title(ax(1),"1N");  
%cum
ax_idx = 1; 
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(cumulativeIN1,2),'.-','Color',colorIN,'MarkerSize',msize+5)
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(cumulativeOUT1,2),'.-','Color',colorOUT,'MarkerSize',msize+5)
%rate
ax_idx = 5; 
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(rateIN1,2),'.-','Color',colorIN,'MarkerSize',msize+5)
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(rateOUT1,2),'.-','Color',colorOUT,'MarkerSize',msize+5)

mask = data_multidic.nf==5&strcmp(data_multidic.spddxl,"lowspd"); 
cumulativeIN2 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
cumulativeOUT2 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
rateIN2 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
rateOUT2 = NaN(max(data_multidic.Nframe)-2,length(find(mask==1))); 
ii = 0; 
for itrial = 1:Ntrial
    subjectNum = find(strcmp(data_multidic.subject(itrial),subjectlist));
    if data_multidic.nf(itrial)==5&&strcmp(data_multidic.spddxl(itrial),"lowspd")
        ii = ii+1; 
        ax_idx = 3;
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JmaxIN,'y')
        %plot(1:data_multidic.Nframe(itrial),results{itrial}.JminOUT,'g')
%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform1+"minOUT"),'--','Color',colorOUT)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform1+"min")(1:end-2,1)*1e2,'.-','Color',[1 1 1].*0.5,'MarkerSize',msize)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform2+"max")(1:end-2,2)*1e2,'.-','Color',[1 1 1].*0.5,'MarkerSize',msize)
        cumulativeIN2(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform1+"min")(1:end-2,1)*1e2;
        cumulativeOUT2(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform2+"max")(1:end-2,2)*1e2;

%         plot(ax(ax_idx),1:data_multidic.Nframe(itrial),results_cum{itrial}.(deform2+"maxIN"),'-','Color',colorIN)
        xlim(ax(ax_idx),[1,data_multidic.Nframe(itrial)]); 
        ylim(ax(ax_idx),[-0.1,0.1]*1e2);

        ax_idx = 7;
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform1+"ratemin")(1:end-2,1)*ratefactor*1e2,'.--','Color',[1 1 1].*0.5,'MarkerSize',msize)
        plot(ax(ax_idx),1:data_multidic.Nframe(itrial)-2,results{itrial}.(deform2+"ratemax")(1:end-2,2)*ratefactor*1e2,'.--','Color',[1 1 1].*0.5,'MarkerSize',msize)
        
        rateIN2(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform1+"ratemin")(1:end-2,1)*ratefactor*1e2;
        rateOUT2(1:data_multidic.Nframe(itrial)-2,ii) = results{itrial}.(deform2+"ratemax")(1:end-2,2)*ratefactor*1e2;

        xlim(ax(ax_idx),[1,data_multidic.Nframe(itrial)]); 
        ylim(ax(ax_idx),[-0.01*ratefactor*1e2,0.01*ratefactor*1e2]);

    end
end

%title 
title(ax(3),"5N");  
%cum
ax_idx = 3; 
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(cumulativeIN2,2),'.-','Color',colorIN,'MarkerSize',msize+5)
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(cumulativeOUT2,2),'.-','Color',colorOUT,'MarkerSize',msize+5)

%rate
ax_idx = 7; 
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(rateIN2,2),'.-','Color',colorIN,'MarkerSize',msize+5)
plot(ax(ax_idx),1:max(data_multidic.Nframe)-2,mean(rateOUT2,2),'.-','Color',colorOUT,'MarkerSize',msize+5)

idx1 = data_multidic.Nframe(data_multidic.nf==1&strcmp(data_multidic.spddxl,"lowspd"))-2;
idx2 = data_multidic.Nframe(data_multidic.nf==5&strcmp(data_multidic.spddxl,"lowspd"))-2;
cIN1 = diag(cumulativeIN1(idx1,:));
cIN2 = diag(cumulativeIN2(idx2,:));
cOUT1 = diag(cumulativeOUT1(idx1,:));
cOUT2 = diag(cumulativeOUT2(idx2,:));

rIN1 = diag(rateIN1(idx1,:));
rIN2 = diag(rateIN2(idx2,:));
rOUT1 = diag(rateOUT1(idx1,:));
rOUT2 = diag(rateOUT2(idx2,:));

h = bar(ax(4),[0+0.25,0+0.75,1.5+0.25,1.5+0.75],[mean(cIN1), mean(cIN2),...
                             mean(cOUT1),mean(cOUT2)],'FaceColor',[1 1 1].*0.5); 
h.ShowBaseLine='off'; 
ax(2).YAxis.Visible = 'on';
bar(ax(8),[0+0.25,0+0.75,1.5+0.25,1.5+0.75],[mean(rIN1), mean(rIN2),...
                             mean(rOUT1),mean(rOUT2)],'FaceColor',[1 1 1].*0.5); 
ax(4).XTick = [0+0.25,0+0.75,1.5+0.25,1.5+0.75]; 
ax(4).XTickLabel = ["1N","5N","1N","5N"]; 
ax(4).XLim = [-0.25,1.5+0.75];
xtickangle(ax(4),60); 
ax(8).XTick = [0+0.25,0+0.75,1.5+0.25,1.5+0.75]; 
ax(8).XTickLabel = ["1N","5N","1N","5N"]; 
ax(8).XLim = [-0.25,1.5+0.75];
xtickangle(ax(8),60)
                         
linkaxes([ax(1:4)],'y'); 
linkaxes([ax(5:8)],'y'); 

set(ax,'FontSize',15); 
set(ax,'Color','None'); 
set(ax(1),'Clipping','off');
text(ax(1),-6,mean(get(ax(1),'YLim')),"Total (%)",'Rotation',90,...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'FontSize',18);

ax(1).XTick = []; 
ax(1).XAxis.Visible = 'off'; 
ax(3).YTickLabel = []; 
ax(3).XAxis.Visible = 'off'; 
ax(3).XTick= [];  
ax(4).YTickLabel = []; 
ax(5).XAxis.Visible = 'off'; 
ax(6).XAxis.Visible = 'off'; 
ax(8).YTickLabel = []; 


title(ax(4),'Final'); 
set(ax(4),'Clipping','off');

set(ax,'Clipping','off');

text(ax(5),-6,mean(get(ax(5),'YLim')),"Rate (%/s)",'Rotation',90,...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'FontSize',18);
ylim6 = get(ax(6),'YLim');


text(ax(4),mean([0+0.25,0+0.75,1.5+0.25,1.5+0.75]),10,'Final',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top',...
    'FontSize',18);

text(ax(6),+6,ylim6(1)-23,"# Frame",'Rotation',0,...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'FontSize',18);
xlabel(ax(6),'# Frame'); 
ax(6).YTickLabel = []; 
ylim(ax(8),[-0.025*ratefactor*1e2,0.01*ratefactor*1e2]);


%X-Axis 
y0 = -10.25;
plot(ax(4-1),[30 (30-12.5)],[y0 y0],'k-','linew',2);
text(ax(4-1),mean([30 (30-12.5)]),y0,"500ms",...
        'FontSize',15,...
        'Rotation',0,...
        'FontWeight','bold',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top');
    
y0 = -61;
plot(ax(end),[30 (30-12.5)],[y0 y0],'k-','linew',2);
text(ax(end-1),mean([30 (30-12.5)]),y0-1,"500ms",...
        'FontSize',15,...
        'Rotation',0,...
        'FontWeight','bold',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top');

%legend
qw = cell(4,1); 
qw{1} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],colorIN);
qw{2} = fill(ax(ax_idx),[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],colorOUT);
qw{3} = plot(ax(ax_idx),NaN,NaN,'-','Color','k');
qw{4} = plot(ax(ax_idx),NaN,NaN,'--','Color','k');
hlegend = legend(ax(3),[qw{:}],{'COMPR.IN','STRET. OUT'});%,'min','max'
set(hlegend, 'NumColumns', 2,...
    'Location','best',...
    'Box','off',...
    'FontSize',15);
box(hlegend,'off'); 
%%
T = results(results.frame ~=0,:);
T = grpstats(T,...
    {'subject','nf'},...
    {'mean'},...
    'DataVars',{'Epc1maxIN','Epc2maxIN','JmaxIN','ShearmaxIN',...
                'Epc1maxOUT','Epc2maxOUT','JmaxOUT','ShearmaxOUT',...
                'Epc1minIN','Epc2minIN','JminIN','ShearminIN',...
                'Epc1minOUT','Epc2minOUT','JminOUT','ShearminOUT'},...
    'VarNames',{'subject','nf','GrpCnt',...
            'Epc1maxIN','Epc2maxIN','JmaxIN','ShearmaxIN',...
            'Epc1maxOUT','Epc2maxOUT','JmaxOUT','ShearmaxOUT',...
            'Epc1minIN','Epc2minIN','JminIN','ShearminIN',...
            'Epc1minOUT','Epc2minOUT','JminOUT','ShearminOUT'}...
); 
figure; 
%subplot(1,2,1); 
hold on; 
subjectlist = unique(T.subject); 

xshift = 0; 
bar([0 1]+xshift,[abs(mean(T.Epc1minIN(T.nf==1))),...
    abs(mean(T.Epc1minIN(T.nf==5)))],'FaceColor',[1 1 1].*0.5); 
for ii = 1:Nsubject
    plot((T.nf(T.subject==subjectlist(ii))==5)+xshift,...
        abs(T.Epc1minIN(T.subject==subjectlist(ii))),'b-'); 
end
%boxplot(abs(T.Epc1minIN),T.nf); 
scatter((T.nf==5)+xshift,abs(T.Epc1minIN),200,'b.'); 
%ylim([0 0.15]);

%subplot(1,2,2); hold on; 
xshift = 3; 
bar([0 1]+xshift,[abs(mean(T.Epc2maxOUT(T.nf==1))),...
    abs(mean(T.Epc2maxOUT(T.nf==5)))],'FaceColor',[1 1 1].*0.5); 
for ii = 1:Nsubject
    plot((T.nf(T.subject==subjectlist(ii))==5)+xshift,...
        abs(T.Epc2maxOUT(T.subject==subjectlist(ii))),'r-'); 
end
scatter((T.nf==5)+xshift,abs(T.Epc2maxOUT),200,'r.');
ylim([0 0.15]); 

%boxplot(abs(T.Epc2maxOUT),T.nf);

%%
z_threshold = 1.1; % height to which we consider to be in contact with the surface 
deform_threshold = 0.5; 
ratedeform_threshold = 0.05; % arbitrary threshold of max deformation between 2 frames 
corrcoef_threshold = 1.5; % arbitrary threshold of max deformation between 2 frames 

ptile = 10; 
Ntrial = length(data_multidic);
Nframe = length(data_multidic(1).Points3D);
Jinmin_all = zeros(Nframe,Ntrial);  
Joutmax_all = zeros(Nframe,Ntrial);  
sum_contact_all = zeros(Nframe,Ntrial);  
for jj = 1:Ntrial
    [def1_in, def1_out, contact_area] = findDeformContact(data_multidic(jj),...
        'deform','J',...
        'z_threshold',z_threshold,...
        'deform_threshold',deform_threshold,...
        'ratedeform_threshold',ratedeform_threshold,...
        'corrcoef_threshold',corrcoef_threshold);
    Nfeat = length(data_multidic(jj).Deform.J{1});
    Jinmean = zeros(Nframe,1); 
    Jinmin = zeros(Nframe,1); 
    Joutmax = zeros(Nframe,1); 
    Joutmean = zeros(Nframe,1); 
    sum_contact = zeros(Nframe,1); 
    for ii = 1:Nframe
        Jinmean(ii) = mean(def1_in{ii},'omitnan'); 
        Jinmin(ii) = prctile(def1_in{ii},ptile); 
        Joutmax(ii) = prctile(def1_out{ii},100-ptile); 
        Joutmean(ii) = mean(def1_out{ii},'omitnan'); 
        sum_contact(ii) = sum(contact_area{ii},'omitnan'); 
    end
    Jinmin_all(:,jj) = Jinmin; 
    Joutmax_all(:,jj) = Joutmax; 
    sum_contact_all(:,jj) = sum_contact; 
end


function varargout = findDeformContact(data,varargin)

% Parse input arguments
p = inputParser;
p.addRequired('data');
p.addParameter('z_threshold',0,@isnumeric);
p.addParameter('deform_threshold',10,@isnumeric);
p.addParameter('ratedeform_threshold',10,@isnumeric);
p.addParameter('corrcoef_threshold',0,@isnumeric);


p.parse(data, varargin{:});
z_threshold = p.Results.z_threshold;
ratedeform_threshold = p.Results.ratedeform_threshold;
deform_threshold = p.Results.deform_threshold;
corrcoef_threshold = p.Results.corrcoef_threshold;

%
Nframe = length(trial.Points3D);

for ii = 1:Nframe
    def = data.dic.Deform(ii,:);
    Fcnow=data.dic.FaceCentroids{ii};
    if ii~=1
        defprevious =  trial.Deform.(deform){ii-1};
        mask1 = abs(def-defprevious)>ratedeform_threshold;
        mask2 = abs(def-1)>deform_threshold;
        def(mask1|mask2,:)=NaN;
    end
    corrcoefmask = trial.FaceCorrComb{ii} < corrcoef_threshold;  
    areamask = trial.FaceCentroids{ii}(:,3) < z_threshold;
    
    %fit plane with point below threshold
    threshold = -0.9; 
    idx_IN = Fcnow(:,3) < threshold; 
    DM = [Pnow(idx_IN,1),Pnow(idx_IN,2),ones(length(Pnow(idx_IN,1)),1)]; % Design Matrix
    param = DM\Pnow(idx_IN,3);                   % Estimate Parameters
    B = param; 
    planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.1);

    %find the boundary 
    plane_id = planefit(Pnow(:,1),Pnow(:,2)) > Pnow(:,3); 
    CAnow = findboundary(Pnow(plane_id,:));
    contactarea{ii} = trial.Deform.Area{ii}(areamask&corrcoefmask,:); 
    def_in{ii} =  def(areamask&corrcoefmask);
    def_out{ii} =  def(~areamask&corrcoefmask);
end
varargout{1} = mask_in; 
end
%%
%refine threshold from last map 
%     threshold = -1; 
%     idx_IN = Fcentroidlast(:,3) < threshold; 
%     DM = [Fcentroidlast(idx_IN,1),Fcentroidlast(idx_IN,2),ones(length(Fcentroidlast(idx_IN,1)),1)]; % Design Matrix
%     param = DM\Fcentroidlast(idx_IN,3);                   % Estimate Parameters
%     B = param; 
%     planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.13);
%     plane_id = planefit(Fcentroidlast(:,1),Fcentroidlast(:,2)) > Fcentroidlast(:,3); 
%     threshold = mean(Fcentroidlast(plane_id,3));
%     
%     %fit plane with point below threshold
%     IN_threshold = 300;
%     idx_IN = Fcentroidnow(:,3) < threshold+0.5; 
%     length(find(idx_IN==1))
%     DM = [Fcentroidnow(idx_IN,1),Fcentroidnow(idx_IN,2),ones(length(Fcentroidnow(idx_IN,1)),1)]; % Design Matrix
%     param = DM\Fcentroidnow(idx_IN,3);                   % Estimate Parameters
%     B = param; 
%     planefit = @(x,y) B(1)*x + B(2)*y + (B(3)+0.13);
% 
%     %find the boundary 
%     plane_id = Fcentroidnow(:,3) < planefit(Fcentroidnow(:,1),Fcentroidnow(:,2)); 
%     if length(find(idx_IN==1))>IN_threshold
%         CAnow = findboundary(Fcentroidnow(plane_id,:));
%     else
%         CAnow = [NaN,NaN,NaN];
%     end