%    complexFreq2afc.m
%
%         by: joshw, june 2021
%
%    purpose: 2afc task comparing images with specific spatial frequency components. 
%             Images are inverse fourier transforms of distributions generated in frequency space.
%             The frequency space images are weighted compilations of gaussian
%             distributed spatial frequency components. 
%             
%             
  

function myscreen = compFreq(varargin)
getArgs(varargin);

% initilaize the screen
myscreen = initScreen('hplp'); mglClearScreen; task{1}.waitForBacktick = 1;

% set task parameters
task{1}.segmin = [2 .5 1 inf]; task{1}.segmax = [2 .5 1 inf]; task{1}.getResponse = [0 0 0 1]; %length of segments and response segment
task{1}.numTrials = 500;
task{1}.parameter.frequencyComp1 = [15]
task{1}.parameter.frequencyComp2 = [30]
task{1}.parameter.weight = [.2 .35 .5 .65 .8]
task{1}.random = 1

% create a multidimensional array of comparison images
for i = 1:50
task{1}.compare(:,:,i) = generateFrequencyImage((i+20)/2,(i+20)/2,.5,myscreen);
end

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task,myscreen)

if task.thistrial.thisseg == 1
    task.thistrial.complexImage = generateFrequencyImage(task.thistrial.frequencyComp1,task.thistrial.frequencyComp2,task.thistrial.weight,myscreen);
    task.thistrial.estimate = 10;
    mglClearScreen(0);myscreen.flushMode = 1;
    %task = jumpSegment(task)
end




if task.thistrial.thisseg == 2
        tex = mglCreateTexture(task.thistrial.complexImage*255);
        mglBltTexture(tex,[0 0]);
        mglFlush
elseif task.thistrial.thisseg == 3
    mglClearScreen(0);mglFlush;
elseif task.thistrial.thisseg == 4
    myscreen.flushMode = 0;
    task.thistrial.gotResponse = 0
    task.thistrial.estimate = 10;task.thistrial.current = 500;
end


mglFixationCross;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
if task.thistrial.thisseg == 4 
    
    tex = mglCreateTexture(task.compare(:,:,task.thistrial.estimate)*255);
    mglBltTexture(tex,[0 0]);
    
    if ~(task.thistrial.estimate == task.thistrial.current);mglFlush;end; %only update if buffer is different
    task.thistrial.current = task.thistrial.estimate;
    
    [task estimateDone] = setEstimate(myscreen,task); %set new comparison for next flush
    
    if estimateDone
        task.thistrial.gotResponse = 1;
        % save trial parameters in a convenient place for analysis (personal preference) IF we get response
        task = jumpSegment(task);
        disp(sprintf('(estimation) Estimate: %0.2f',task.thistrial.estimate))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%
%       estimate       %
%%%%%%%%%%%%%%%%%%%%%%%%
function [task estimateDone] = setEstimate(myscreen,task)

% get scroll
scrollEvents = mglListener('getAllScrollEvents');
if ~isempty(scrollEvents)
  % get the sum of the vertical and horizontal scrolls
  verticalScroll = sum(scrollEvents.scrollVertical);
  horizontalScroll = sum(scrollEvents.scrollHorizontal);
  % set estimate
  task.thistrial.estimate = task.thistrial.estimate+verticalScroll;;
else
  horizontalScroll = 0;
end


% check bounds
if task.thistrial.estimate > 50,task.thistrial.estimate = 50;end;if task.thistrial.estimate < 2,task.thistrial.estimate = 2;end;

% if mouse button down (or horizontal scroll is non-zero) then we are done setting estimate
mouse = mglGetMouse;
keyEvent = mglGetKeyEvent;

if ~isequal(mouse.buttons,0) || ~isempty(keyEvent)
  estimateDone = 1;
else
  estimateDone = 0;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)
if task.thistrial.whichButton == 1
   task = jumpSegment(task);
end
if task.thistrial.whichButton == 2
   task = jumpSegment(task);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% generateFrequencyImages %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img = generateFrequencyImage(frequency1,frequency2,weight,myscreen)

spatialFreq = frequency1; dev = (spatialFreq/3);
spatialFreq2 = frequency2; dev2 = (spatialFreq2/3);
weight = weight;

stimSize = 1000;

dTheta = 0.01; dTheta2 = 0.01;
nSamples = 1/dTheta; nSamples2 = 1/dTheta2;

%% generate stimulus in the frequency domain
% generating a semi-circle

freqStimComplex = zeros(stimSize,stimSize);

for radius = 1:max(spatialFreq,spatialFreq2)*2.5

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
