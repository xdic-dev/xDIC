
function [frameRangeStart, frameRangeEnd] = calculateFrameRanges(phase, NbrFrLoading, NbrFrSlide, NbrFrRelax)
% Calculate the frame ranges based on the specified phase

switch phase 
    case 'loading' 
        frameRangeStart = 1;
        NbrFrame = NbrFrLoading; 
        frameRangeEnd = frameRangeStart+NbrFrame-1;
    case 'slide1'
        frameRangeStart = NbrFrLoading+1;
        NbrFrame = NbrFrSlide; 
        frameRangeEnd = frameRangeStart+NbrFrame-1;
    case 'relax1'
        frameRangeStart = NbrFrLoading+NbrFrSlide+1;
        NbrFrame = NbrFrRelax; 
        frameRangeEnd = frameRangeStart+NbrFrame-1;
    case 'slide2'
        frameRangeStart = NbrFrLoading+NbrFrSlide+NbrFrRelax+1;
        NbrFrame = NbrFrSlide; 
        frameRangeEnd = frameRangeStart+NbrFrame-1;
    case 'relax2'
        frameRangeStart = NbrFrLoading+2*NbrFrSlide+NbrFrRelax+1;
        NbrFrame = NbrFrRelax; 
        frameRangeEnd = frameRangeStart+NbrFrame-1;
    case 'all'
        frameRangeStart = 1; 
        frameRangeEnd = NbrFrLoading+2*NbrFrSlide+2*NbrFrRelax; 
    otherwise
        error('Invalid phase. Please provide a valid phase type.');
end
end