%% probabilistic population coding with poisson like neural variability


% inputs %
tuningStd=10 %Gaussian tuning curves at the moment
numNeurons=200
boundary=150 %boundary of the charts from 0 to X
numStims=900 %Number of points in PLD

stimulus=50
gain=200
stimulus2=100
gain2=200

r = zeros([1,numNeurons]);
stimVal = zeros([1,numStims]);

integrate = 1; %% if you want to run simulations with multiple inputs

%%%% simulate neural activity with poisson variability
for i = 1:numNeurons
    avgRate = gain*normpdf(stimulus,i*(boundary/numNeurons),tuningStd);
    r(i)=poissrnd(avgRate);
end
figure(1);subplot(1,2,1);scatter(((boundary/numNeurons):(boundary/numNeurons):boundary),r);title('Population Response');
xlabel('Preferred orientation');ylabel('Spikes/second');

for s = boundary/numStims:boundary/numStims:boundary
    p = 1;
    for ri = 1:length(r);
        p = p*poisspdf(r(ri),gain*normpdf(s,ri*(boundary/numNeurons),tuningStd));
    end
    stimVal(round(s*(numStims/boundary))) = p;
end
stimVal = stimVal*(1/sum(stimVal));
subplot(1,2,2);scatter(boundary/numStims:boundary/numStims:length(stimVal)/(numStims/boundary),stimVal);title('Stimulus Probability');
title('Perceptual likelihood distribution'),xlabel('Location'),ylabel('Likelihood');

cumulative = zeros([1,length(stimVal)]);

for i = 2:length(cumulative)
    cumulative(i) = cumulative(i-1)+stimVal(i);
end
PLDfit = fitCumulativeGaussian([boundary/numStims:boundary/numStims:length(cumulative)/(numStims/boundary)],cumulative)
L1 = sprintf('Stimulus: %d // Mean: %.2f // Std: %.2f',stimulus,PLDfit.mean,PLDfit.std)
legend(L1)
xlim([30,115])



if integrate %%%% if you want to integrate other stimuli %%%%%
   
r2 = zeros([1,numNeurons]);
stimVal2 = zeros([1,numStims]);

    for i = 1:numNeurons
        avgRate2 = gain2*normpdf(stimulus2,i*(boundary/numNeurons),tuningStd);
        r2(i)=poissrnd(avgRate2);
    end
    figure(2);subplot(1,2,1);scatter(((boundary/numNeurons):(boundary/numNeurons):boundary),r2);title('Population Response');
    xlabel('Preferred orientation');ylabel('Spikes/second');

    for s = boundary/numStims:boundary/numStims:boundary
        p = 1;
        for ri = 1:length(r2);
            p = p*poisspdf(r2(ri),gain2*normpdf(s,ri*(boundary/numNeurons),tuningStd));
        end
        stimVal2(round(s*(numStims/boundary))) = p;
    end
    stimVal2 = stimVal2*(1/sum(stimVal2));
    figure(2);subplot(1,2,2);scatter(boundary/numStims:boundary/numStims:length(stimVal2)/(numStims/boundary),stimVal2);title('Stimulus Probability');
    title('Perceptual likelihood distribution'),xlabel('Location'),ylabel('Likelihood');

    cumulative2 = zeros([1,length(stimVal2)]);

    for i = 2:length(cumulative2)
        cumulative2(i) = cumulative2(i-1)+stimVal2(i);
    end
    PLDfit2 = fitCumulativeGaussian([boundary/numStims:boundary/numStims:length(cumulative2)/(numStims/boundary)],cumulative2)
    L2 = sprintf('Stimulus: %d // Mean: %.2f // Std: %.2f',stimulus2,PLDfit2.mean,PLDfit2.std)
    legend(L2)
    xlim([30,115])
   
%%%% combine %%%%

%r3 neural activity
r3 = r+r2;
figure(3);subplot(1,2,1);scatter(((boundary/numNeurons):(boundary/numNeurons):boundary),r3);title('Population Response');
xlabel('Preferred orientation');ylabel('Spikes/second');

% likelihood
stimVal3 = zeros([1,numStims]);
for s = boundary/numStims:boundary/numStims:boundary
    p = 1;
    for ri = 1:length(r3);
        p = p*poisspdf(r3(ri),(gain+gain2)*normpdf(s,ri*(boundary/numNeurons),tuningStd));
    end
    stimVal3(round(s*(numStims/boundary))) = p;
end
stimVal3 = stimVal3*(1/sum(stimVal3));
figure(3);subplot(1,2,2);scatter(boundary/numStims:boundary/numStims:length(stimVal3)/(numStims/boundary),stimVal3);title('Stimulus Probability');
hold on;scatter(boundary/numStims:boundary/numStims:length(stimVal2)/(numStims/boundary),stimVal2,'red');scatter(boundary/numStims:boundary/numStims:length(stimVal)/(numStims/boundary),stimVal,'green');
title('Perceptual likelihood distribution'),xlabel('Location'),ylabel('Likelihood');

cumulative3 = zeros([1,length(stimVal3)]);

 for i = 2:length(cumulative3)
     cumulative3(i) = cumulative3(i-1)+stimVal3(i);
 end
    PLDfit3 = fitCumulativeGaussian([boundary/numStims:boundary/numStims:length(cumulative3)/(numStims/boundary)],cumulative3)
   
    %%% statistics %%%
    std1 = PLDfit.std;std2=PLDfit2.std;std3=PLDfit3.std;
    m1 = PLDfit.mean;m2=PLDfit2.mean;m3=PLDfit3.mean;
    pstd = sqrt(((std1^2)*(std2^2))/(std1^2+std2^2))
    pm = m1*((std2^2)/(std1^2+std2^2))+m2*((std1^2)/(std1^2+std2^2))

    L3 = sprintf('Mean: %.2f // Std: %.2f // Predicted Mean: %.2f // Predicted Std: %.2f', PLDfit3.mean,PLDfit3.std,pm,pstd)
    legend(L3,L1,L2)
    xlim([30,115])
end