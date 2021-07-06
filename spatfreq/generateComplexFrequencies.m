% generateStim.m - generate half of a stimulus in the frequency domain and
% inverse fourier transform it back into the spatial domain

clear;plot=1;


spatialFreq = 30; dev = (spatialFreq/3);
spatialFreq2 = 10; dev2 = (spatialFreq2/3);
weight = .9

stimSize = 1000;

dTheta = 0.001; dTheta2 = 0.001;
nSamples = 1/dTheta; nSamples2 = 1/dTheta2;

%% generate stimulus in the frequency domain
% generating a semi-circle

freqStimComplex = zeros(stimSize,stimSize);

for radius = 1:100%max(spatialFreq,spatialFreq2)*2.5

for i = 0:nSamples-1
    xPlot(i+1) = round(radius*cos(pi*dTheta*i))+round(stimSize/2); 
    yPlot(i+1) = -1*(round(radius*sin(pi*dTheta*i)))+round(stimSize/2); 
end

for i = 0:nSamples2-1
    xxPlot(i+1) = round(radius*cos(pi*dTheta2*i))+round(stimSize/2);
    yyPlot(i+1) = -1*(round(radius*sin(pi*dTheta2*i)))+round(stimSize/2);
end

theta = 0:2*pi/length(xPlot):2*pi-1/length(xPlot);
theta = theta(randperm(length(theta)));

theta2 = 0:2*pi/length(xxPlot):2*pi-1/length(xxPlot);
theta2 = theta2(randperm(length(theta2)));

for j = 1:length(xPlot)
    % theta(j) = 2*pi*rand;
    freqStimComplex(xPlot(j),yPlot(j)) = weight*normpdf(radius,spatialFreq,dev)*(cos(theta(j)) + sin(theta(j))*1i);
end

for j = 1:length(xxPlot)
    freqStimComplex(xxPlot(j),yyPlot(j)) = freqStimComplex(xxPlot(j),yyPlot(j)) + (1-weight)*normpdf(radius,spatialFreq2,dev2)*(cos(theta2(j)) + sin(theta2(j))*1i);
end
end

%% getting the other half of the frequency domain and inverse transforming it to the spatial domain

halfStimFreq = getHalfFourierA(freqStimComplex);
[img,imgf] = reconstructFromHalfFourierA(halfStimFreq);

img=img-min(img(:)); % shift data such that the smallest element of A is 0
img=img/max(img(:)); % normalize the shifted data to 1 










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
colormap(gray);
colorbar
title('Stimulus (Spatial)')
axis off;
axis square
end