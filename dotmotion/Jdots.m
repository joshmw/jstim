% Jdots.m
%
%       $Id$
%       usage: Jdots (temporary)
%       by: josh wilson
%       date: March 2021
%       copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%       
%
%    Purpose: Dot motion. Can parameterize up to 7 motion components per trial & change noise type and motion transparency.
%
%    Things you can change: 
%
%       Dot Speed: set 'stimulus.speed', in degrees/second.
%
%       Dot Density: set 'stimulus.density', in dots/degree.
%
%       Motion Direction and Coherence: This program needs *7 inputs*. You can either define 7 different direction
%       and coherence parameters (i.e. 'task.parameter.coherence1 = [.2 .3]', 'task.parameter.coherence2 = [.15]', etc) and
%       the program will build the array input for you or you can parameterize and pass in your own pre-built 7-term arrays
%       to stimulus.dots.setDir and stimulus.dots.setCoherence. The Nth terms of the direction (in circular degrees) and
%       coherence (percentage of dots that will move in that direction each frame update) arrays correspond. 
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
myscreen = initScreen('hplp');
mglClearScreen
task{1}.waitForBacktick = 1;

global stimulus

% task parameters
task{1}.segmin = [5 1 7 1];
task{1}.segmax = [5 1 7 1];
task{1}.getResponse = [0 0 1 0];
task{1}.numTrials = 100;
task{1}.randVars.calculated.confidence = nan;
% stimulus values are constant throughout experiment
stimulus.movshonNoise = 1; %set either movshon noise OR random walk to 1
stimulus.randomWalk = 0;
stimulus.transparent = 1;
stimulus.Ecc = 0;
stimulus.width = 10;
stimulus.speed = 4;
stimulus.contrast = .8;
stimulus.density = 2.5;
% component 
task{1}.parameter.dir = [45];
task{1}.parameter.dir2 = [55];
task{1}.parameter.dir3 = [35];
task{1}.parameter.dir4 = [60];
task{1}.parameter.dir5 = [30];
task{1}.parameter.dir6 = [65];
task{1}.parameter.dir7 = [25];
task{1}.parameter.coherence = [.15];
task{1}.parameter.coherence2 = [.3];
task{1}.parameter.coherence3 = [.3];
task{1}.parameter.coherence4 = [.1];
task{1}.parameter.coherence5 = [.1];
task{1}.parameter.coherence6 = [.025];
task{1}.parameter.coherence7 = [.025];

% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback);
end

% init the stimulus
if (stimulus.movshonNoise & stimulus.randomWalk) | ~(stimulus.movshonNoise | stimulus.randomWalk); disp('!!!!!!! Too many or too few noise types set !!!!!!!');keyboard;end;
if max(task{1}.parameter.coherence)+max(task{1}.parameter.coherence2)+max(task{1}.parameter.coherence3)+max(task{1}.parameter.coherence4)+max(task{1}.parameter.coherence5+max(task{1}.parameter.coherence6)+max(task{1}.parameter.coherence7)) > 1; disp(' !!!!!! Cannot have cumulative coherence > 1 !!!!!!');keyboard;end;
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);


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
    stimulus.dots = stimulus.dots.setCoherence(stimulus.dots,[task.thistrial.coherence,task.thistrial.coherence2,task.thistrial.coherence3,task.thistrial.coherence4,task.thistrial.coherence5,task.thistrial.coherence6,task.thistrial.coherence7]);
    
    % set direction
    stimulus.dots = stimulus.dots.setDir(stimulus.dots,[task.thistrial.dir,task.thistrial.dir2,task.thistrial.dir3,task.thistrial.dir4,task.thistrial.dir5,task.thistrial.dir6,task.thistrial.dir7]);
    
    
    % send initdots whether you want a transparents or nontransparent motion
    stimulus.dots = stimulus.dots.setTransparency(stimulus.dots,stimulus.transparent,rand(1,stimulus.width*stimulus.width*stimulus.density))
    % I've coded this so that it the dots.transparent term == 1, the dotvals will be based on the transparentDotvals array
    % instead of randomizing each trial, effectively keeping each dot in the same coherence for each frame update instead of randomizing
    % which coherence it belongs. the last term (width*width*density) = dots.n in the initdots script.
    
    % set the speed
    %stimulus.speed(task.trialnum) = stimulus.speed;
    %stimulus.dots = stimulus.dots.setSpeed(stimulus.dots,stimulus.speed(task.trialnum));
    
    % set the fixation color
    stimulus.fixColor = stimulus.stimulatingFixColor;
    
elseif (task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4)

    stimulus.dots = stimulus.dots.setCoherence(stimulus.dots,[0 0 0 0 0 0 0]);
    
    % set the fixation color to indicate response period
    stimulus.fixColor = stimulus.responseFixColor;
    

elseif task.thistrial.thisseg == 3
  % set starting confidence
  task.thistrial.confidence = 0;
  task.thistrial.gotResponse = 0;
  scrollEvents = mglListener('getAllScrollEvents');
  mglListener('getAllMouseEvents');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

% clear screen to gray
mglClearScreen(stimulus.backgroundColor);

if task.thistrial.thisseg == 1
    
    % update the dots
    stimulus.dots = stimulus.dots.update(stimulus.dots);
    
    % draw the dots
    stimulus.dots = stimulus.dots.draw(stimulus.dots);
    
elseif task.thistrial.thisseg == 2 | task.thistrial.thisseg == 4
    
    % update the dots
    stimulus.dots = stimulus.dots.update(stimulus.dots);
    
    % draw the dots
    stimulus.dots = stimulus.dots.draw(stimulus.dots);
    
end

if task.thistrial.thisseg == 3
  % set the confidence
  [task.thistrial.confidence confidenceDone] = setConfidence(task.thistrial.confidence, stimulus);
  if confidenceDone
    task.thistrial.gotResponse = 1;
    % save trial parameters in a convenient place for analysis (personal preference) IF we get response
    task = jumpSegment(task);
    disp(sprintf('(estimation) Estimate: %0.2f',task.thistrial.confidence))
  end
end

% draw fixation cross
mglFixationCross(1,1,stimulus.fixColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = getResponseCallback(task,myscreen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dots stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen,task)

width = strcat('width=',num2str(stimulus.width));
contrast = strcat('contrast=',num2str(stimulus.contrast));
speed = strcat('speed=',num2str(stimulus.speed));
Ecc = strcat('xCenter=',num2str(stimulus.Ecc));
dotDense = strcat('density=',num2str(stimulus.density));
movshon = strcat('movshonNoise=',num2str(stimulus.movshonNoise));
randwalk = strcat('randomWalk=',num2str(stimulus.randomWalk));

% init the dot patches
stimulus.dots = initJdots('framesPerSecond',myscreen.framesPerSecond,width,contrast,speed,Ecc,dotDense,movshon,randwalk);

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
function [confidence confidenceDone] = setConfidence(confidence, stimulus)

% get scroll
scrollEvents = mglListener('getAllScrollEvents');
if ~isempty(scrollEvents)
  % get the sum of the vertical and horizontal scrolls
  verticalScroll = sum(scrollEvents.scrollVertical);
  horizontalScroll = -sum(scrollEvents.scrollHorizontal);
  % set confidence
  confidence = confidence+verticalScroll;
else
  horizontalScroll = 0;
end

% check bounds
if confidence > 360,confidence = 0;end
if confidence < 0,confidence = 360;end

% draw the confidence
drawConfidence(confidence);

% if mouse button down (or horizontal scroll is non-zero) then we are done setting confidence
mouse = mglGetMouse;
keyEvent = mglGetKeyEvent;

if ~isequal(mouse.buttons,0) || ~isempty(keyEvent)
  confidenceDone = 1;
else
  confidenceDone = 0;
end

function drawConfidence(confidenceLevel)
guessValue = sprintf('%.0f',confidenceLevel*100);
%mglTextDraw(guessValue,[0 0]);
mglLines2(0,0,7*cos(deg2rad(confidence)),7*sin(deg2rad(confidence)),3,[0 0 0])



