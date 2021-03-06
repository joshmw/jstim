function [fullImage noiseImage gaussianImage] = generateBackground
myscreen = initScreen;

%%% set parameters %%%
width = [.5 2.5]
weight = .75 %% weight of the first width in forming the amplitude. weight of second is automatically 1-weight.
SNR = 0
maxSNR = 1
contrast = 1
  
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

for iWidth = 1:length(width)
    
% get the gaussian full size of the screen

if iWidth == 1
[Background(iWidth).gaussian Background(iWidth).x Background(iWidth).y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, width(iWidth), width(iWidth));
else
[Background(iWidth).gaussian Background(iWidth).x Background(iWidth).y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, width(iWidth), width(iWidth), 5,2);
end


% now make sure that dimensions are odd numbers of pixels
oddWidth = 2*floor(myscreen.screenWidth/2)+1;
oddHeight = 2*floor(myscreen.screenHeight/2)+1;

% resize everything to odd
Background(iWidth).gaussian = Background(iWidth).gaussian(1:oddHeight,1:oddWidth);
Background(iWidth).x = Background(iWidth).x(1:oddHeight,1:oddWidth);
Background(iWidth).y = Background(iWidth).y(1:oddHeight,1:oddWidth);
  
% get the fourier transform
Background(iWidth).gaussianTransform = getHalfFourier(Background(iWidth).gaussian);
  
% pull out magnitude and dc for averaging
mag = Background(iWidth).gaussianTransform.mag;
dc = Background(iWidth).gaussianTransform.dc;

% make average transform
Background(iWidth).averageGaussianTransform = Background(iWidth).gaussianTransform;
Background(iWidth).averageGaussianTransform.dc = dc;
Background(iWidth).averageGaussianTransform.mag = mag;

end


%%%%% Reconstruct the image %%%%%%
% max noise and signal
noiseMax = 1 / (maxSNR + 1);
sigMax = SNR * noiseMax;

% randomize phase and reconstruct
Background(3).averageGaussianTransform.phase = (rand(1,(Background(1).averageGaussianTransform.n+Background(2).averageGaussianTransform.n)/2)*2*pi - pi);
Background(3).averageGaussianTransform.dc = (Background(1).averageGaussianTransform.dc + Background(2).averageGaussianTransform.dc)/2;

% here we weight the amplitude inputs of the component images
Background(3).averageGaussianTransform.mag = (weight*Background(1).averageGaussianTransform.mag + (1-weight)*Background(2).averageGaussianTransform.mag)/2;

Background(3).averageGaussianTransform.originalDims = Background(1).averageGaussianTransform.originalDims;

Background(1).averageGaussianTransform.phase = (rand(1,Background(1).averageGaussianTransform.n)*2*pi - pi);
Background(1).im = reconstructFromHalfFourier(Background(1).averageGaussianTransform);

Background(2).averageGaussianTransform.phase = (rand(1,Background(2).averageGaussianTransform.n)*2*pi - pi);
Background(2).im = reconstructFromHalfFourier(Background(2).averageGaussianTransform);

Background(3).im = reconstructFromHalfFourier(Background(3).averageGaussianTransform);

%%% gaussian %%%
Background(3).gaussianTransform = Background(1).gaussianTransform;
Background(3).gaussianTransform.mag = (weight*Background(1).gaussianTransform.mag+(1-weight)*Background(2).gaussianTransform.mag)/2;
Background(3).gaussian = reconstructFromHalfFourier(Background(3).gaussianTransform);



for i = 1:3
% scale from 0 to noise max
maxIm = max(Background(i).im(:));
minIm = min(Background(i).im(:));
Background(i).im(:,:) = noiseMax * (Background(i).im - minIm) / (maxIm-minIm);
end


%plot
for i = 1:3
figure(i);imshow(Background(i).im);
end

mglClose
k=2
keyboard










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
