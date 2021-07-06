function image = imgsynth
global stimulus
  iWidth=1
  myscreen = initScreen
  maxIndex = 255;
  stimulus.colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];
  stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
  stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
  stimulus.colors.nGaussianColors = maxIndex+1-stimulus.colors.nReservedColors;
  stimulus.colors.minGaussianIndex = maxIndex+1 - stimulus.colors.nGaussianColors;
  stimulus.colors.midGaussianIndex = stimulus.colors.minGaussianIndex + floor(stimulus.colors.nGaussianColors/2);
  stimulus.colors.maxGaussianIndex = maxIndex;
  stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nGaussianColors-1);
  stimulus.colors.reservedColors = [1 1 1; 0.3 0.3 0.3; 0 1 0;1 0 0; 0 1 1];
  stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
  
  setGammaTableForMaxContrast(.5);
  contrastIndex = getContrastIndex(.5,1);
  
  
  
  
  % get the gaussian full size of the screen
  [G.gaussian{iWidth} G.x G.y] = mglMakeGaussian(myscreen.imageWidth, myscreen.imageHeight, 1, 1);

  % now make sure that dimensions are odd numbers of pixels
  oddWidth = 2*floor(myscreen.screenWidth/2)+1;
  oddHeight = 2*floor(myscreen.screenHeight/2)+1;

  % resize everything to odd
  G.gaussian{iWidth} = G.gaussian{iWidth}(1:oddHeight,1:oddWidth);
  G.x = G.x(1:oddHeight,1:oddWidth);
  G.y = G.y(1:oddHeight,1:oddWidth);
  
  % get the fourier transform
  G.gaussianTransform{iWidth} = getHalfFourier(G.gaussian{iWidth});
  
  % pull out magnitude and dc for averaging
  mag(iWidth,:) = G.gaussianTransform{iWidth}.mag;
  dc(iWidth) = G.gaussianTransform{iWidth}.dc;

% make average transform
G.averageGaussianTransform = G.gaussianTransform{1};
G.averageGaussianTransform.dc = dc;
G.averageGaussianTransform.mag = mag(3:3:end);





G.averageGaussianTransform.phase = (rand(1,length(G.averageGaussianTransform.mag))*2*pi - pi);
  im = reconstructFromHalfFourier(G.averageGaussianTransform);

  % scale from 0 to noise max
  maxIm = max(im(:));
  minIm = min(im(:));
  G.im(1,:,:) = G.noiseMax * (im - minIm) / (maxIm-minIm);

  % make into texture
  G.texture(iBackground) = mglCreateTexture(round(stimulus.colors.gaussRange*squeeze(G.im(1,:,:)) + stimulus.colors.minGaussianIndex));
  
  % update disppercent
  disppercent(iBackground/numBackgrounds);
  
  
  
  
  
  
  
  
  
  function setGammaTableForMaxContrast(maxContrast)

global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

% set the reserved colors
gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% set the gamma table
if maxContrast > 0
  % create the rest of the gamma table
%   cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
  cmin = 0;
  cmax = maxContrast;
  luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nGaussianColors-1)):cmax;

  % replace NaN in gamma tables with zero
  stimulus.linearizedGammaTable.redTable(isnan(stimulus.linearizedGammaTable.redTable)) = 0;
  stimulus.linearizedGammaTable.greenTable(isnan(stimulus.linearizedGammaTable.greenTable)) = 0;
  stimulus.linearizedGammaTable.blueTable(isnan(stimulus.linearizedGammaTable.blueTable)) = 0;

  % now get the linearized range
  redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  % add these values to the table
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
  % if we are asked for 0 contrast then simply set all the values to BLACK
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0,'linear');
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0,'linear');
  gammaTable((stimulus.colors.minGaussianIndex:stimulus.colors.maxGaussianIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% keep the gamma table
stimulus.gammaTable = gammaTable;

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;



function contrastIndex = getContrastIndex(desiredContrast,verbose)

if nargin < 2,verbose = 0;end

global stimulus;
if desiredContrast < 0, desiredContrast = 0;end

% now find closest matching contrast we can display with this gamma table
contrastIndex = min(round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.colors.nDisplayContrasts);

% display the desired and actual contrast values if verbose is set
if verbose
  actualContrast = stimulus.currentMaxContrast*(contrastIndex/stimulus.colors.nDisplayContrasts);
  disp(sprintf('(getContrastIndex) Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,actualContrast,desiredContrast-actualContrast));
end

% out of range check
if round(stimulus.colors.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.colors.nDisplayContrasts
 disp(sprintf('(getContrastIndex) Desired contrast (%0.9f) out of range max contrast : %0.9f',desiredContrast,stimulus.currentMaxContrast));
 keyboard
end

% 1 based indexes (0th index is gray, nDisplayContrasts+1 is full contrast)
contrastIndex = contrastIndex+1;
