%%%%%%%%%%%%%%%%%%%%%%%%
%    getHalfFourier    %
%%%%%%%%%%%%%%%%%%%%%%%%
function d = getHalfFourier(im)

% make sure there are an odd number of pixels
if mod(size(im,1),2) == 0
    im = im(1:end-1,:);
end
if mod(size(im,2),2) == 0
    im = im(:,1:end-1);
end

% take fourier transform of image
%imf = fft2(im);
imf = im;

% get input dimensions
d.originalDims = size(im);

% get one half of fourier image
% all of the x-dim, first half of the y-dim (e.g. imfh(1:401,1:201);)
%imfh = fftshift(imf);
imfh = imf(1:d.originalDims(1),1:ceil(d.originalDims(2)/2));


% extract dc from half fourier image
d.dc = imfh(ceil(d.originalDims(1)/2),end);
halfFourier = imfh(1:(prod(size(imfh))-ceil(d.originalDims(1)/2)));

d.mag = abs(halfFourier);
d.phase = angle(halfFourier);
d.n = length(d.phase);
