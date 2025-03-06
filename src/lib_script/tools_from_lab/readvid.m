%% readvid
% 
% reads images from a compressed h264 mp4 video
%
% inputs:
%    file        path to the mp4 video 
%    npic        [optional] number of pictures to read (default=all)
%    jump        [optional] number of pictures to skip before reading (default=0)
%    every       [optional] read 1 frame every x frames (default=1)
%
% output:
%    imarray     array of pictures, dim = (m,n,npic)



function imarray=readvid(file,npic,jump,every)

v=VideoReader(file);

if(nargin<2),   npic=round(v.Duration*v.FrameRate);end
if(nargin<3),   jump=0;end
if(nargin<4),   every=1;end

s=[v.Height v.Width npic];

% handle jumps
%v.CurrentTime=(jump+1)/v.FrameRate;

imarray=zeros(s,'uint8');
% for kk=1:npic
%   for jj=1:every
%     iloc=readFrame(v);
%   end
%   imarray(:,:,kk)=iloc(:,:,1);
% end
for kk=1:npic
  iloc=read(v,jump+(kk-1)*every+1);
  imarray(:,:,kk)=iloc(:,:,1);
end

end

