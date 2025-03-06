function out = bandpassfft(in,r1,r2,debug)
%BANDPASSFFT Band-pass filter on images in Fourier domain.
%Gaussian band pass filter in Fourier domain for images.
%
% Syntax: out = bandpassfft(in,r1,r2,debug)
%
% Inputs:
%   in        image or image array (of class double or uint8)
%   r1        small radius of the band pass filter (high-pass threshold)
%   r2        large radius of the band pass filter (low-pass threshold)
%   debug     (bool) debug mode if 1
%
% Outputs:
%   out       image or image array after filtering (of class double)
%
% Example:
%   close all; clear all; clc;
%   im = imread('fingerprint.tif');
%   im_filtered = bandpassfft(im,20,200,1);
%   figure; imshow(im_filtered,[]);
%
% Note: to understand how to use fftshift and ifftshift:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/285244
%
% A band-pass filter is the combination of a low-pass filter and a high-
% pass filter, where the cut-off frequency of the low-pass is higher than
% that of the high pass. There is always a trade-off between blurring and
% noise: low-pass reduces noise but accentuates blurring, high-pass reduces
% blurring but accentuates noise.
% Source: http://scribd.com/doc/51981950/Frequency-Domain-Bandpass-Filtering-for-Image-Processing
%
% Note: to visualize the filtered images properly, use imshow(im,[]).

% Process input
if(nargin < 4)
    debug = 0;
end

if(r1 >= r2)
    error('r1 must be smaller than r2!');
end

% Create gaussian filter
siz = size(in);
c = ceil(siz/2);
xgv = (1:siz(2)) - c(2);
ygv = (1:siz(1)) - c(1);
[x,y] = meshgrid(xgv,ygv);
r = sqrt(x.^2+y.^2);
if(r1 ~= 0)
    filter = (exp(-r.^2/(2*r2^2)).*(1-exp(-r.^2/(2*r1^2))));
else
    filter = exp(-r.^2/(2*r2^2));
end

% If ~debug, avoid creating multiple huge matrices for memory usage
if debug
    % Compute 2D DFT of input image
    F = fft2(ifftshift(in));
    % Apply filter
    I = bsxfun(@times,F,ifftshift(filter));
    % Compute 2D inverse DFT
    out = real(fftshift(ifft2(I)));
else
    % All in one command without ifftshift on the input/output.
    % Same result, because, not using fftshift / ifftshift only introduces
    % errors on the phase of the signal, which doesn't matter with images.
    out = real(ifft2(fft2(in).*fftshift(filter)));
end

% (DEBUG) Display in, image DFT (both magnitude and phase), filter and out
if debug
    % Input and output images
    figure;
    subplot(2,1,1);
    imshow(in,[]); title('Input image');
    subplot(2,1,2);
    imshow(out,[]); title('Output image');

    % Filter
    figure; surf(imresize(filter,0.1)); title('Filter');

    % Image DFT before and after filtering (with both magnitude and phase)
    figure;
    subplot(2,2,1);
    M = visualize_magnitude_fft(F);
    imagesc(M); colormap(gray);
    title('Magnitude spectrum BEFORE filtering');

    subplot(2,2,3);
    A = visualize_angle_fft(F);
    imagesc(A); colormap(gray);
    title('Phase spectrum BEFORE filtering');

    subplot(2,2,2);
    M = visualize_magnitude_fft(I);
    imagesc(M); colormap(gray);
    title('Magnitude spectrum AFTER filtering');

    subplot(2,2,4);
    A = visualize_angle_fft(I);
    imagesc(A); colormap(gray);
    title('Phase spectrum AFTER filtering');
end

    %----------------------------------------------------------------------
    function M = visualize_magnitude_fft(F)
        % Source: http://stackoverflow.com/questions/13549186/how-to-plot-a-2d-fft-in-matlab
        M = fftshift(F); % Center DFT
        M = abs(M);      % Get the magnitude
        M = log(M+1);    % Use log, for perceptual scaling, and +1 since log(0) is undefined
        M = 100*M;       % Arbitrary scaling factor to improve image display
        M = mat2gray(M); % Use mat2gray to scale the image between 0 and 1
    end

    %----------------------------------------------------------------------
    function A = visualize_angle_fft(F)
        A = fftshift(F); % Center DFT
        A = angle(A);    % Get the phase
        A = mat2gray(A); % Use mat2gray to scale the image between 0 and 1
    end

end
