%%% simulate a 2afc vs estimation perceptual likelihood distribution %%%
PLFstd = 3; %std of the cue
CLFstd = 3; %std of the 2afc comparison distribution
ELFstd = 3; %std of the estimation comparison distribution
start = -15;finish = 15;

%%%%% make underlying perceptual distribution %%%%%
x = [start:0.1:finish];
%PLF = normpdf(x,0,PLFstd); kl = sum(PLF);PLF = PLF*(1/sum(PLF)); %normal
PLF = normpdf(x,-5,PLFstd)+normpdf(x,5,PLFstd); k = sum(PLF);PLF = PLF*(1/sum(PLF)); %bimodal
 

%%%%% simulate 2afc draws %%%%%
Offsets = [-15 -12 -9 -7 -5 -3 -1 0 1 3 5 7 9 12 15]; % probe offsets
numRepeats = 40; % number trials per offset
Rightward = ones(1,length(Offsets));

for offset = 1:length(Offsets)
for repeat = 2:numRepeats  
   
% make comparison distribution for 2afc (probe cue) - can try making different ones
CLF = normpdf(x,Offsets(offset),CLFstd); k = sum(CLF);CLF = CLF*(1/sum(CLF)); %unimodal
%CLF = normpdf(x,Offsets(offset)-5,CLFstd)+normpdf(x,Offsets(offset)+5,CLFstd); k = sum(CLF);CLF = CLF*(1/sum(CLF)); %bimodal

% draw from each
PLFPercept = rand(1); CLFPercept = rand(1); PLFComp = 0; CLFComp = 0;
for j = 1:length(PLF)
    if PLFComp > PLFPercept
        PLFX = x(j);break;
    else PLFComp = PLFComp+PLF(j);end;
end
for j = 1:length(CLF)
    if CLFComp > CLFPercept
        CLFX = x(j);break;
    else CLFComp = CLFComp+CLF(j);end;
end

% compare
if PLFX < CLFX
    Rightward(offset) = Rightward(offset)+1;
end
end
end
Rightward = Rightward/numRepeats


%%%%%% simulate estimations %%%%%%
Estimates = []
numEsts = numRepeats*5

for p = 1:numEsts
   
% Draw from PLF
PLFPercept = rand(1); PLFComp = 0; ELFPercept = rand(1); ELFComp = 0;
for j = 1:length(PLF)
    if PLFComp > PLFPercept
        PLFX = x(j);break;
    else PLFComp = PLFComp+PLF(j);end;
end
 
% Draw from estimate distribution
ELF = normpdf(x,PLFX,ELFstd); k = sum(ELF);ELF = ELF*(1/sum(ELF)); %make and normalize to 1 - can change type of distribution?
for j = 1:length(ELF)
    if ELFComp > ELFPercept
        ELFX = x(j); Estimates = [Estimates ELFX]; break;
    else ELFComp = ELFComp+ELF(j);end;
end
end

 






%%%%% figures %%%%%%
q=2
if q==2
    figure(1)
    scatter(x,PLF*numEsts*length(PLF)/20)
    title('Theoretical percetual likelihood distribution')
    xlabel('Arbitrary Units')
    ylabel('Likelihood')
   
   
    figure(2)
    hold off
    scatter(Offsets,Rightward,'blue')
    z=fitCumulativeGaussian(Offsets,Rightward)
    hold on
    scatter([start:.01:finish],normcdf([start:.01:finish],z.mean,z.std),.1,'blue')
    scatter([start:.01:finish],normpdf([start:.01:finish],z.mean,z.std)*10,.1,'black')
    xlabel('Probe offset (a.u.)')
    ylabel('Likehood percieved >0')
    H = sprintf('2afc Prediction: mu = %.2f // std = %.2f // R2 = %.2f',z.mean,z.std,z.r2)
    title(H)
   
     
 
    figure(3)
    hist(Estimates,20)
    xlabel('Estimate')
    ylabel('Frequency')
    H = sprintf('Estimated value: mu = %.2f // std = %.2f',mean(Estimates),std(Estimates))
    title(H)
end