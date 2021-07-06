function [fullImage noiseImage gaussianImage] = generateBackground
myscreen = initScreen;

%%% set parameters %%%
width = [6 4 2]
SNR = 1
maxSNR = 1
contrast = 10
centerWidth = 16
  
%%% set colors, gamma table, contrast index %%%
global Colors; global stimulus;
xPos = 0; yPos = 0;
Colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];
Colors.nReservedColors = size(Colors.reservedColors,1);
maxIndex = 255;
Colors.nGaussianColors = maxIndex+1-Colors.nReservedColors;
Colors.minGaussianIndex = maxIndex+1 - Colors.nGaussianColors;
Colors.maxGaussianIndex = maxIndex;
stimulus.linearizedGammaTable = mglGetGammaTable;
Colors.nDisplayContrasts = floor(Colors.nGaussianColors-1);
setGammaTableForMaxContrast(contrast);
contrastIndex = getContrastIndex(contrast,1);
Colors.gaussRange = contrastIndex-1;
 
for iFreq = 1:length(width)

% get the gaussian full size of the screen
[Background.gaussian Background.x Background.y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, width(iFreq), width(iFreq));

% now make sure that dimensions are odd numbers of pixels
oddWidth = 2*floor(myscreen.screenWidth/2)+1;
oddHeight = 2*floor(myscreen.screenHeight/2)+1;

% resize everything to odd
Background.gaussian = Background.gaussian(1:oddHeight,1:oddWidth);
Background.x = Background.x(1:oddHeight,1:oddWidth);
Background.y = Background.y(1:oddHeight,1:oddWidth);
  
% get the fourier transform
Background.gaussianTransform = getHalfFourier(Background.gaussian);
  
% pull out magnitude and dc for averaging
mag = Background.gaussianTransform.mag;
dc = Background.gaussianTransform.dc;

% make average transform
Background.averageGaussianTransform = Background.gaussianTransform;
Background.averageGaussianTransform.dc = dc;
Background.averageGaussianTransform.mag = mag;

% max noise and signal
Background.noiseMax = 1 / (maxSNR + 1);
Background.sigMax = SNR * Background.noiseMax;

% randomize phase and reconstruct
Background.averageGaussianTransform.phase = (rand(1,Background.averageGaussianTransform.n)*2*pi - pi);
im = reconstructFromHalfFourier(Background.averageGaussianTransform);

% scale from 0 to noise max
maxIm = max(im(:));
minIm = min(im(:));
Background.im(:,:) = Background.noiseMax * (im - minIm) / (maxIm-minIm);

% make into texture
BackTexture = mglCreateTexture(round(Colors.gaussRange*squeeze(Background.im(:,:)) + Colors.minGaussianIndex));

% Images
fullImage = squeeze(Background.im(:,:)) + Background.sigMax * exp(-((((Background.x-xPos).^2) + (Background.y+yPos).^2))/(2*(width(iFreq)^2)));
noiseImage = Background.im;
gaussianImage = Background.gaussian;

% scale and set on texture // unnecessary right now
stimulus.stimTexture = mglCreateTexture(round(Colors.gaussRange*im + Colors.minGaussianIndex));

fullFT = fft2(fullImage(2:end,2:end)-mean(fullImage(:)));
fullCenter = fftshift(fullFT);
fullCenter = fullCenter((-centerWidth/2 : centerWidth/2) +1+size(fullCenter,1)/2,(-centerWidth/2 : centerWidth/2)+1+size(fullCenter,2)/2);

noiseFT = fft2(noiseImage(2:end,2:end)-mean(noiseImage(:)));
noiseCenter = fftshift(noiseFT);
noiseCenter = noiseCenter((-centerWidth/2 : centerWidth/2) +1+size(noiseCenter,1)/2,(-centerWidth/2 : centerWidth/2)+1+size(noiseCenter,2)/2);

gaussianFT = fft2(gaussianImage(2:end,2:end)-mean(gaussianImage(:)));
gaussianCenter = fftshift(gaussianFT);
gaussianCenter = gaussianCenter((-centerWidth/2 : centerWidth/2) +1+size(gaussianCenter,1)/2,(-centerWidth/2 : centerWidth/2)+1+size(gaussianCenter,2)/2);

figure(iFreq);
subplot(2,3,1); imagesc(abs(gaussianCenter));title(sprintf('gaussian: width=%d', width(iFreq)));
subplot(2,3,2); imagesc(abs(noiseCenter)); title('noise');
subplot(2,3,3); imagesc(abs(fullCenter)); title('gaussian + noise');

figure(iFreq);
subplot(2,3,4); imagesc(angle(gaussianCenter));title(sprintf('gaussian: width=%d', width(iFreq)));
subplot(2,3,5); imagesc(angle(noiseCenter)); title('noise');
subplot(2,3,6); imagesc(angle(fullCenter)); title('gaussian + noise');


end
mglClose
k=2


%%%%%%%%%%%%%%%%%%%%%%
%% getContrastIndex %%
%%%%%%%%%%%%%%%%%%%%%%
function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

global stimulus; global Colors;
if desiredContrast < 0, desiredContrast = 0;end

% now find closest matching contrast we can display with this gamma table
contrastIndex = min(round(Colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),Colors.nDisplayContrasts);

% display the desired and actual contrast values if verbose is set
if verbose
  actualContrast = stimulus.currentMaxContrast*(contrastIndex/Colors.nDisplayContrasts);
  disp(sprintf('(getContrastIndex) Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,actualContrast,desiredContrast-actualContrast));
end

% out of range check
if round(Colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>Colors.nDisplayContrasts
 disp(sprintf('(getContrastIndex) Desired contrast (%0.9f) out of range max contrast : %0.9f',desiredContrast,stimulus.currentMaxContrast));
 keyboard
end

% 1 based indexes (0th index is gray, nDisplayContrasts+1 is full contrast)
contrastIndex = contrastIndex+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setGammaTableForMaxContrast %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTableForMaxContrast(maxContrast)

global Colors; global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

% set the reserved colors
gammaTable(1:size(Colors.reservedColors,1),1:size(Colors.reservedColors,2))=Colors.reservedColors;

% set the gamma table
if maxContrast > 0
  % create the rest of the gamma table
  cmin = 0;
  cmax = maxContrast;
  luminanceVals = cmin:((cmax-cmin)/(Colors.nGaussianColors-1)):cmax;

  % replace NaN in gamma tables with zero
  stimulus.linearizedGammaTable.redTable(isnan(stimulus.linearizedGammaTable.redTable)) = 0;
  stimulus.linearizedGammaTable.greenTable(isnan(stimulus.linearizedGammaTable.greenTable)) = 0;
  stimulus.linearizedGammaTable.blueTable(isnan(stimulus.linearizedGammaTable.blueTable)) = 0;

  % now get the linearized range
  redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  % add these values to the table
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
  % if we are asked for 0 contrast then simply set all the values to BLACK
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0,'linear');
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0,'linear');
  gammaTable((Colors.minGaussianIndex:Colors.maxGaussianIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% keep the gamma table
stimulus.gammaTable = gammaTable;

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;
