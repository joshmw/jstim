% generateStim.m - generate half of a stimulus in the frequency domain and
% inverse fourier transform it back into the spatial domain

clear

plot=1
save=0

spatialFreq = 10;
scatter = 5 % controls hexoginality of the image - if set equal to spatialFreq, will be hexagons. lower = more noisy.
stimSize = 100;
dTheta = 0.01;
nSamples = 1/dTheta;
numimages= 1

%% generate stimulus in the frequency domain
for imagenum = 1:numimages
% generating a semi-circle

for i = 0:nSamples-1
    
    xPlot(i+1) = round(spatialFreq*cos(pi*dTheta*i))+round(stimSize/2); % +randi([-round(spatialFreq/scatter) round(spatialFreq/scatter)]);
    yPlot(i+1) = -1*(round(spatialFreq*sin(pi*dTheta*i)))+round(stimSize/2); % +randi([-round(spatialFreq/scatter) round(spatialFreq/scatter)]);
    
end

freqStimComplex = zeros(stimSize,stimSize);

theta = 0:2*pi/length(xPlot):2*pi-1/length(xPlot);
theta = theta(randperm(length(theta)));

for j = 1:length(xPlot)
    % theta(j) = 2*pi*rand;
    freqStimComplex(xPlot(j),yPlot(j)) = (cos(theta(j)) + sin(theta(j))*1i);
end

%% getting the other half of the frequency domain and inverse transforming it to the spatial domain


halfStimFreq = getHalfFourierA(freqStimComplex);
[img,imgf] = reconstructFromHalfFourierA(halfStimFreq);




if save
img=img-min(img(:)); % shift data such that the smallest element of A is 0
img=img/max(img(:)); % normalize the shifted data to 1 

name = ['image' num2str(imagenum) '.png']
imwrite(img,name)
end
end








%% plotting
if plot == 1
fplotRange = spatialFreq*3;
freqRange = ceil(stimSize/2);
figure(1)
imagesc([-freqRange:freqRange],[-freqRange,freqRange],abs(freqStimComplex))
colormap(gray);
colorbar
title('Stimulus (Frequency)')
xlabel('Spatial Frequency X (cycles/image)');
ylabel('Spatial Frequency Y (cycles/image)');
xlim([-fplotRange fplotRange])
ylim([-fplotRange fplotRange])
axis square

% with complex stimulus in the frequency domain
figure(2)
imagesc([-freqRange:freqRange],[-freqRange,freqRange],abs(imgf))
colormap(gray);
colorbar
title('Stimulus (Frequency)')
xlabel('Spatial Frequency X (cycles/image)');
ylabel('Spatial Frequency Y (cycles/image)');
xlim([-fplotRange fplotRange])
ylim([-fplotRange fplotRange])
axis square

figure(3)
imagesc(img);
colormap gray;
colorbar
title('Stimulus (Spatial)')
axis off;
axis square
end