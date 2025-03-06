function DIC3DPPresults = filter_measures(DIC3DPPresults,freqFilt,freqAcq,plotOptList)

    for ii = 1:length(plotOptList)%1:length(optStruct.plotoptdeform)
        % J adjustment
        if strcmp(plotOptList{ii},'J')
            data = DIC3DPPresults.Deform.(plotOptList{ii});
            Nframe = length(data); 
            for kk = 1:Nframe
                DIC3DPPresults.Deform.(plotOptList{ii}){kk} = DIC3DPPresults.Deform.(plotOptList{ii}){kk}-1;
            end
        end
        if strcmp(plotOptList{ii},'DispX')
            data = DIC3DPPresults.Disp.DispVec;
            data_filtered = filter3Ddeform_time(data,...
                'freqFilt',freqFilt,...
                'freqAcq',freqAcq,...
                'deformCol',1);
            DIC3DPPresults.Disp.DispVec = data_filtered;
        elseif strcmp(plotOptList{ii},'DispY')
            data = DIC3DPPresults.Disp.DispVec;
            data_filtered = filter3Ddeform_time(data,...
                'freqFilt',freqFilt,...
                'freqAcq',freqAcq,...
                'deformCol',2);
            DIC3DPPresults.Disp.DispVec = data_filtered;
        elseif strcmp(plotOptList{ii},'DispZ')
            data = DIC3DPPresults.Disp.DispVec;
            data_filtered = filter3Ddeform_time(data,...
                'freqFilt',freqFilt,...
                'freqAcq',freqAcq,...
                'deformCol',3);
            DIC3DPPresults.Disp.DispVec = data_filtered;
        elseif strcmp(plotOptList{ii},'DispMgn')
            data = DIC3DPPresults.Disp.DispMgn;
            data_filtered = filter3Ddeform_time(data,...
                'freqFilt',freqFilt,...
                'freqAcq',freqAcq);
            DIC3DPPresults.Disp.DispMgn = data_filtered;
        elseif strcmp(plotOptList{ii},'FaceCorrComb')
            data = DIC3DPPresults.FaceCorrComb;
            data_filtered = filter3Ddeform_time(data,...
                'freqFilt',freqFilt,...
                'freqAcq',freqAcq);
            DIC3DPPresults.FaceCorrComb = data_filtered;
        else
            data = DIC3DPPresults.Deform.(plotOptList{ii});
            data_filtered = filter3Ddeform_time(data,...
                'freqFilt',freqFilt,...
                'freqAcq',freqAcq);
            DIC3DPPresults.Deform.(plotOptList{ii}) = data_filtered;
        end
    end
end