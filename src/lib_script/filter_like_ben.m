function [imbdp_norm,y] = filter_like_ben(imarray,varargin)

    p = inputParser;
    p.addParameter('paramfilt',[20, 200]);
    p.addParameter('mask',[]);
    p.addParameter('gsbound',[]);
    p.parse(varargin{:});
    paramfilt=p.Results.paramfilt;
    mask=p.Results.mask;
    y=p.Results.gsbound;

%     imarray = cam_mask1; 
    siz = size(imarray); 
    % Filter image array with bandpassfft
    r1=paramfilt(1); r2=paramfilt(2);
    imbdp = zeros(siz);
%     y_all = zeros(1,2); 
    mask = logical(mask);
%         y_all(1) = y_all(1)+y(1); 
%         y_all(2) = y_all(2)+y(2); 

    for ii = 1:siz(3)
        imbdp(:,:,ii) = bandpassfft(imarray(:,:,ii),r1,r2);
        if ii == 1 && isempty(y)
            imbdp_ii = imbdp(:,:,ii); 
            y = prctile(imbdp_ii(mask),[5 95],'all');
        end
%         mask_ii = ncorr1.current_save(ii).roi.mask;
        % Remove extreme values from images
    end
%     y_all = y_all/siz(3);
    imbdp_norm = ((imbdp-y(1))/(y(2)-y(1))); % useful information is from 0 to 1
    imbdp_norm = satur(satur(imbdp_norm,'method','low','level',0),...
        'method','high','level',1);
end