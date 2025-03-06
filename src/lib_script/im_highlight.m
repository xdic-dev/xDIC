
function im = im_highlight(im,mask,varargin)
p = inputParser;
p.addRequired('im');
p.addRequired('mask');
p.addParameter('factor',2,@isnumeric);

p.parse(im,mask,varargin{:});
factor = p.Results.factor;
im(~mask) = im(~mask)*1/factor; 
end