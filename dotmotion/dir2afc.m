%    Jdot2Patch.m
%
%       $Id$
%       usage: Jdot2Patch
%       by: josh wilson
%       date: April 2021
%       copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%       
%
%    Purpose: Dot motion with estimate. Can parameterize up to 7 motion components per patch.
%             Can change noise type and motion transparency, patch location, individual patch components.
%             Can make estimates of speed or direction via scroll wheel in response segment.
%
%    Things you can change: 
%
%       Estimate you're making: stimulus.estimateSpeed = 1 or stimulus.estimateDirection = 1.   
%       Estimates are made by updating an estimate dot patch (stimulus.dots.Response) via scroll wheel.
%
%       Dot Density: set 'stimulus.density', in dots/degree.
%
%       Dot Speed: set 'task{1}.parameter.speed', in degrees/second.
%
%       Motion Direction and Coherence: This program needs *7 inputs*. You can either define 7 different direction
%       and coherence parameters (i.e. 'task.parameter.coherence1 = [.2 .3]' ... 'task.parameter.coherence7 = [.15]', etc) and
%       the program will build the array input for you or you can parameterize and pass in your own pre-built 7-term arrays
%       to stimulus.dots.setDir and stimulus.dots.setCoherence. The Nth terms of the direction (in circular degrees) and
%       coherence (percentage of dots that will move in that direction each frame update) arrays correspond. 
%       Currently set up to just define 1 motion.
%    
%       Motion Transparency: 'stimulus.transparency = 1' or = 0. Which component a dot belongs to can be updated either every
%       frame (nontransparent) or every trial (transparent). Transparent motion looks like multiple gratings sliding
%       over eachother where nontransparent motion does not.
%
%       Noise Type: Can specify movshon or random walk - 'stimulus.movshonNoise = 1' or 'stimulus.randomWalk = 1'.




function myscreen = dotDiscrep(varargin)

% check arguments
getArgs(varargin,'stimulusType=dots');


% initilaize the screen
myscreen = initScreen('hplp'); mglClearScreen; task{1}.waitForBacktick = 1;
global stimulus


% task parameters
task{1}.segmin = [.3 .5 .3 .5 inf 1]; task{1}.segmax = [.3 .5 .3 .5 inf 1]; task{1}.getResponse = [0 0 0 0 1 0];
task{1}.numTrials = 500; task{1}.randVars.calculated.estimate = nan; task{1}.randVars.calculated.dir = nan;


% Set noise and transparency %
stimulus.movshonNoise = 0; %set either movshon noise OR random walk to 1
stimulus.randomWalk = 1;
stimulus.transparent = 0; %nontransparency randomizes which dots belong with each motion component every frame update
stimulus.numStim = 1 %number of stimuli you want to draw
task{1}.random=1;


% Estimate: Set what parameter you want to estimate %
stimulus.estimateSpeed = 0;
stimulus.estimateDirection = 1;


% Set stimulus values %
stimulus.Ecc1 = 8;
stimulus.Ecc2 = -8;
stimulus.width = 8; %stimuli and response
stimulus.speed = 3; %initial response speed
stimulus.contrast = .6; %stimuli and response
stimulus.density = 2; %stimuli and response


% Stimuli motion components
task{1}.parameter.dir1 = [1:360];
task{1}.parameter.coherence1 = [.3];
task{1}.parameter.speed = [3];

task{1}.parameter.dir2 = [140];
task{1}.parameter.coherence2 = [.6];
task{1}.parameter.speed2 = [5];

task{1}.parameter.offset = [-20 -15 -10 -7 -4 -2 .001 2 4 7 10 15 20]

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end

% init the stimulus
if (stimulus.movshonNoise & stimulus.randomWalk) | ~(stimulus.movshonNoise | stimulus.randomWalk); disp('!!!!!!! Too many or too few noise types set !!!!!!!');keyboard;end;
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

% response variables
task{1}.response.correct = nan
task{1}.response.diff = nan

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
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus
if task.thistrial.thisseg == 1
        
    % set the coherence
    stimulus.dots = stimulus.dots.setCoherence(stimulus.dots,[task.thistrial.coherence1,0,0,0,0,0,0]);
    stimulus.dots2 = stimulus.dots2.setCoherence(stimulus.dots2,[task.thistrial.coherence2,0,0,0,0,0,0]);
    
    % set direction
    stimulus.dots = stimulus.dots.setDir(stimulus.dots,[task.thistrial.dir1,0,0,0,0,0,0]);
    stimulus.dots2 = stimulus.dots2.setDir(stimulus.dots2,[task.thistrial.dir2,0,0,0,0,0,0]);
    stimulus.dotsResponse = stimulus.dotsResponse.setDir(stimulus.dotsResponse,[task.thistrial.dir2,0,0,0,0,0,0]);
    
    % send initdots whether you want a transparents or nontransparent motion
    stimulus.dots = stimulus.dots.setTransparency(stimulus.dots,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density));
    stimulus.dots2 = stimulus.dots2.setTransparency(stimulus.dots2,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density));
    stimulus.dotsResponse = stimulus.dotsResponse.setTransparency(stimulus.dotsResponse,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density));
    % I've coded this so that it the dots.transparent term == 1, the dotvals will be based on the transparentDotvals array
    % instead of randomizing each trial, effectively keeping each dot in the same coherence for each frame update instead of randomizing
    % which coherence it belongs. the last term (width*width*density) = dots.n in the initdots script.
    
    % set the speed
    stimulus.dots = stimulus.dots.setSpeed(stimulus.dots,task.thistrial.speed);
    stimulus.dots2 = stimulus.dots.setSpeed(stimulus.dots2,task.thistrial.speed2);
    
    % set the fixation color
    stimulus.fixColor = stimulus.stimulatingFixColor;
    
elseif task.thistrial.thisseg == 3
    
    % set the coherence
    stimulus.dots = stimulus.dots.setCoherence(stimulus.dots,[task.thistrial.coherence1,0,0,0,0,0,0]);
    stimulus.dots2 = stimulus.dots2.setCoherence(stimulus.dots2,[task.thistrial.coherence2,0,0,0,0,0,0]);
    
    % set direction
    dir = task.thistrial.dir1+task.thistrial.offset
    if dir > 360; dir = dir-360; end;
    if dir > 0; dir = dir+360; end;
    stimulus.dots = stimulus.dots.setDir(stimulus.dots,[dir,0,0,0,0,0,0]);
    stimulus.dots2 = stimulus.dots2.setDir(stimulus.dots2,[dir,0,0,0,0,0,0]);
    
    stimulus.dotsResponse = stimulus.dotsResponse.setDir(stimulus.dotsResponse,[task.thistrial.dir2,0,0,0,0,0,0]);
    
    % send initdots whether you want a transparents or nontransparent motion
    stimulus.dots = stimulus.dots.setTransparency(stimulus.dots,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density));
    stimulus.dots2 = stimulus.dots2.setTransparency(stimulus.dots2,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density));
    stimulus.dotsResponse = stimulus.dotsResponse.setTransparency(stimulus.dotsResponse,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density));
    % I've coded this so that it the dots.transparent term == 1, the dotvals will be based on the transparentDotvals array
    % instead of randomizing each trial, effectively keeping each dot in the same coherence for each frame update instead of randomizing
    % which coherence it belongs. the last term (width*width*density) = dots.n in the initdots script.
    
    % set the speed
    stimulus.dots = stimulus.dots.setSpeed(stimulus.dots,task.thistrial.speed);
    stimulus.dots2 = stimulus.dots.setSpeed(stimulus.dots2,task.thistrial.speed2);
    
    % set the fixation color
    stimulus.fixColor = stimulus.stimulatingFixColor;
    
    
elseif (task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4 | task.thistrial.thisseg == 6)

    stimulus.dots = stimulus.dots.setCoherence(stimulus.dots,[0 0 0 0 0 0 0]);
    stimulus.dots2 = stimulus.dots2.setCoherence(stimulus.dots2,[0 0 0 0 0 0 0]);
    
    % set the fixation color to indicate response period
    stimulus.fixColor = stimulus.responseFixColor;
    

elseif task.thistrial.thisseg == 5
  % set starting confidence
  task.thistrial.gotResponse = 0;
  scrollEvents = mglListener('getAllScrollEvents');
  mglListener('getAllMouseEvents');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%resp%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1
    
    % update the dots
    stimulus.dots = stimulus.dots.update(stimulus.dots);
    stimulus.dots2 = stimulus.dots2.update(stimulus.dots2);
    
    % draw the dots
    stimulus.dots = stimulus.dots.draw(stimulus.dots);
    if stimulus.numStim>1
    stimulus.dots2 = stimulus.dots2.draw(stimulus.dots2);
    end
    
elseif task.thistrial.thisseg == 3
    % update the dots
    stimulus.dots = stimulus.dots.update(stimulus.dots);
    stimulus.dots2 = stimulus.dots2.update(stimulus.dots2);
    
    % draw the dots
    stimulus.dots = stimulus.dots.draw(stimulus.dots);
    if stimulus.numStim>1
    stimulus.dots2 = stimulus.dots2.draw(stimulus.dots2);
    end
    
elseif task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4 | task.thistrial.thisseg == 6
    
    % update the dots
    stimulus.dots = stimulus.dots.update(stimulus.dots);
    stimulus.dots2 = stimulus.dots2.update(stimulus.dots2);
    
    % draw the dots
    stimulus.dots = stimulus.dots.draw(stimulus.dots);
    if stimulus.numStim>1
    stimulus.dots2 = stimulus.dots2.draw(stimulus.dots2);
    end
end

if task.thistrial.thisseg == 5
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)
if task.thistrial.whichButton == 1
   if task.thistrial.offset > 0; task.response.correct = [task.response.correct 1];end
   if task.thistrial.offset < 0; task.response.correct = [task.response.correct 0];end
   task.response.diff = [task.response.diff task.thistrial.offset]
   task = jumpSegment(task);
end
if task.thistrial.whichButton == 2
   if task.thistrial.offset > 0; task.response.correct = [task.response.correct 0];end
   if task.thistrial.offset < 0; task.response.correct = [task.response.correct 1];end
   task.response.diff = [task.response.diff task.thistrial.offset]
   task = jumpSegment(task);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen,task)

width = strcat('width=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.speed));
Ecc1 = strcat('xCenter=',num2str(stimulus.Ecc1));
Ecc2 = strcat('xCenter=',num2str(stimulus.Ecc2));
dotDense = strcat('density=',num2str(stimulus.density));
movshon = strcat('movshonNoise=',num2str(stimulus.movshonNoise));
randwalk = strcat('randomWalk=',num2str(stimulus.randomWalk));

% init the dot patches
stimulus.dots = initJdots('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,Ecc1,dotDense,movshon,randwalk);
stimulus.dots2 = initJdots('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,Ecc2,dotDense,movshon,randwalk);
stimulus.dotsResponse = initJdots('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,'xCenter=0',dotDense,movshon,randwalk);


% set background color
stimulus.backgroundColor = 0.5;

% and fixation color
stimulus.stimulatingFixColor = [0 0 0];
stimulus.responseFixColor = [1 1 1];
stimulus.correctFixColor = [0 1 0];
stimulus.incorrectFixColor = [1 0 0];
stimulus.neutralFixColor = [0 0 1];


%%%%%%%%%%%%%%%%%%%%%%%%
%       estimate       %
%%%%%%%%%%%%%%%%%%%%%%%%
function [confidence confidenceDone] = setConfidence(confidence,stimulus)

% get scroll
scrollEvents = mglListener('getAllScrollEvents');
if ~isempty(scrollEvents)
  % get the sum of the vertical and horizontal scrolls
  verticalScroll = sum(scrollEvents.scrollVertical);
  horizontalScroll = -sum(scrollEvents.scrollHorizontal);
  % set confidence
  if stimulus.estimateSpeed==1;confidence = confidence+verticalScroll/50;end;
  if stimulus.estimateDirection==1;confidence = confidence+verticalScroll;end;
else
  horizontalScroll = 0;
end


% check bounds
if stimulus.estimateSpeed==1;
    if confidence > 10,confidence = 10;end;if confidence < 0,confidence = 0;end;
elseif stimulus.estimateDirection==1;
    if confidence > 360,confidence = 0;end;if confidence < 0,confidence = 360;end;
end

% if mouse button down (or horizontal scroll is non-zero) then we are done setting confidence
mouse = mglGetMouse;
keyEvent = mglGetKeyEvent;

if ~isequal(mouse.buttons,0) || ~isempty(keyEvent)
  confidenceDone = 1;
else
  confidenceDone = 0;
end









k=1
if k ==2
diff = task{1}.response.diff(2:end); correct = task{1}.response.correct(2:end);
percents = []
for i = 1:length(task{1}.parameter.offset)
numcor = 0;total = 0;
for k = 1:500
if diff(k) == task{1}.parameter.offset(i); numcor = numcor+correct(k);total = total+1;end
end
percents = [percents numcor/total]
end
for i = 1:6
percents(i) = 1-percents(i)
end
figure(100);scatter(task{1}.parameter.offset,percents)
z=fitCumulativeGaussian(task{1}.parameter.offset,percents)
hold on
scatter([-30:.01:30],normcdf([-30:.01:30],z.mean,z.std),.1,'black')
xlabel('Probe offset (degrees)');ylabel('Percentage seen leftward)');
H = sprintf('2afc direction of motion: mu = %.2f // std = %.2f',z.mean,z.std)
title(H)
end
